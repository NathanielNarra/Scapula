 function [EL_final, SS] = Discrete_Ricci_Flow_G0(TR)

% Calculate the Riemannian metric
[RM, Edge_Lengths] = Discrete_Riemannian_Metric(TR);

% Calculate the initial circle metric
[I, Gamma_f, Gamma_v,I_list]= Initial_Circle_Packing_Metric(TR,RM,Edge_Lengths);

[EL_final, SS] = Newton_Method(TR,I,Gamma_f,Gamma_v,I_list);

% Embedding onto a planar domain
% Flat_V = Seed_face(EL_final,TR);
%Flat_V = Planar_Embedding(EL_final,TR,SS.CornerAngles(:,:,end)); % doesnt
%work

 end

function [EW, Gamma_f, Gamma_v, InvD] = Initial_Circle_Packing_Metric(TR,RM,EL)
TRI = TR.ConnectivityList;

% Compute Gamma values for each vertex for each Triangle (Face)
Gamma_f = (repmat(sum(RM,2),[1,3])-2*RM)/2;

% Compute Gamma values for each vertex as a function of its multiple values
% from adjacent triangles
% Gamma_v = struct2array(regionprops(TRI,Gamma_f,'MinIntensity'))';
Gamma_v = ones(size(TR.Points,1),1);
% Compute inversive distance (RECHECK)
temp = Gamma_v(edges(TR));
InvD = (EL.^2 - sum(temp.^2,2))./(2*prod(temp,2));
% fun = @(x) Gamma_v(x.data);
% G = blockproc(TRI,[1 1],fun);
Gamma_f = Gamma_v(TRI);
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]); 
EW = (RM.^2 - Gi.^2 - Gj.^2)./(2*Gi.*Gj);
%EW(EW<0) = 0;
end

function [EL, SS] = Newton_Method(TR,I,~,Gamma_v,I_list)
noV = size(TR.Points,1);
U = log(Gamma_v);
E = edges(TR);
SS = {};
%epsilon = 0.05;
SS.U0 = U;

% U = U - mean(U);    
Gamma_v = exp(U);
Gamma_f = Gamma_v(TR.ConnectivityList);


%**** Optimse below
FE = Faces_constituting_Edges(TR);
%[Edge_indices_in_F, Edge_indices_in_F_mask] = Edge_Instances(TR.ConnectivityList,FE);
Edge_indices_in_F = Edge_Instances(TR,FE);
[V_matrix_edgeid,~,V_Occur,~] = Vertex_Connectivity(TR);
%[Vertex_all_edges,Vertex_all_edges_mask]  = Edges_for_each_Vertex(E,noV);
%****


% Current curvature
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]);
EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2);
denominator = repmat(2*prod(EL,2),[1,3])./EL;
corner_angles =  acos(numerator./denominator);

% Optimised code below
% Theta_sums = sum((corner_angles(V_Occur)).*V_Occur_mask)';
Theta_sums = full(sum(spfun(@(x) corner_angles(x),V_Occur))');
% angle_list = regionprops(TR.ConnectivtyList,Gamma_f,'PixelValues');
% Theta_sums = cellfun(@sum,struct2cell(angle_list))';

K = 2*pi - Theta_sums;
% Mixed Area
% Area_mixed = Mesh_MixedArea(EL,corner_angles,TR.ConnectivityList);
% Current curvature
% K = K./Area_mixed;
% Target Curvature
% K_target = (4*pi)/sum(Area_mixed);
K_target = ones(size(Gamma_v))*(4*pi)/size(TR.Points,1);

Error = max(abs(K_target - K));
iteration_count = 1;
while Error > 0.0000001 & iteration_count<150
    
    % Edge lengths
    Gi = circshift(Gamma_f,[0 -1]);
    Gj = circshift(Gamma_f,[0 -2]);
    EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
    
    % Computing corner angles
    numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2); % INSERT CHECK FOR ERRONEOUS NEGATIVE VALUES. sign of too acute angles/triangles
    denominator = repmat(2*prod(EL,2),[1,3])./EL;
    Q = numerator./denominator;
    if max(Q) >1 | min(Q)<-1
        trimesh(TR), axis equal
        hold on
        ix = find(any(Q>1,2));
        p1 = plot3(TR.Points(TR.ConnectivityList(ix,:),1),TR.Points(TR.ConnectivityList(ix,:),2),TR.Points(TR.ConnectivityList(ix,:),3),'*r');
        hold off
        error('acos out of bounds: probably bad triangle quality')
    end
    corner_angles =  acos(numerator./denominator);
    Area_mixed = Mesh_MixedArea(EL,corner_angles,TR.ConnectivityList);
    Area_sums(iteration_count,1)= sum(Area_mixed);

    % Current curvature
    % Optimised code below
%     Theta_sums = sum((corner_angles(V_Occur)).*V_Occur_mask)';
    Theta_sums = full(sum(spfun(@(x) corner_angles(x),V_Occur))');
%     angle_list = regionprops(TR.ConnectivtyList,Gamma_f,'PixelValues');
%     Theta_sums = cellfun(@sum,struct2cell(angle_list))';
    %
    
    %Area_mixed = Mesh_MixedArea(EL,corner_angles,TR.ConnectivityList);
%     K = (2*pi - Theta_sums)./Area_mixed;
    K = (2*pi - Theta_sums);
    
    % Parameters to detect scale space features
    SS.U(:,iteration_count) = U;
    Gij = Gamma_v(E);
    SS.L(:,iteration_count) = sqrt(sum(Gij.^2,2) + (prod(Gij,2).*I_list));
    SS.K(:,iteration_count) = K;
    SS.CornerAngles(:,:,iteration_count) = corner_angles;
    SS.EL(:,:,iteration_count) = EL;    
    
    %H = Hessian_Matrix(TR,EL,Gamma_f,I,Edge_indices_in_F,Edge_indices_in_F_mask,V_matrix_edgeid,V_matrix_edgeid_mask,V_Occur,V_Occur_mask);
    % Orthocentre and distances therefrom 
    Orth_dist = Orthocentre_Edge_Distance(TR,EL,corner_angles,Gamma_f);
    
    %Edge weights and Hessian Matrix construction.
    W_matrix = Orth_dist./EL;
    W_edge_weights = full(sum(spfun(@(x) W_matrix(x),Edge_indices_in_F))');
    H = spfun(@(x) W_edge_weights(x),V_matrix_edgeid);
    H = diag(sum(H,2)) - H;

    % Preconditioning and Conjugate gradient method
    L = ichol(H);    
    dU = pcg(H,(K_target - K),[],noV,L,L');
    U = U + dU; % all node metrics are updated
    U = U - mean(U);
    
    Gamma_v = exp(U);
    Gamma_f = Gamma_v(TR.ConnectivityList);
    
    % Update target curvature
%     K_target = (4*pi)/sum(Area_mixed);
    
    %Error(iteration_count+1) = max(abs(K_target-K));
    Error = max(abs(K_target-K))
    iteration_count = iteration_count + 1;
        

end


%% Final edge lengths
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]);
EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2);
denominator = repmat(2*prod(EL,2),[1,3])./EL;
corner_angles =  acos(numerator./denominator);
Theta_sums = full(sum(spfun(@(x) corner_angles(x),V_Occur))');
K = (2*pi - Theta_sums);

SS.U(:,iteration_count) = U;
Gij = Gamma_v(E);
SS.L(:,iteration_count) = sqrt(sum(Gij.^2,2) + (prod(Gij,2).*I_list));
SS.K(:,iteration_count) = K;
SS.CornerAngles(:,:,iteration_count) = corner_angles;
SS.EL(:,:,iteration_count) = EL;
end