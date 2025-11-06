function verification_error = verify_validation_orbit(Orbit,validation_terms_one,validation_terms_two,r_ddot)
h = Orbit.h;
h_dot = Orbit.h_dot;

num_time_points = size(h,2);
num_h_modes = size(h,1);

h_ddot_one = zeros(num_h_modes,num_time_points);
h_ddot_two = zeros(num_h_modes,num_time_points);
verification_error_all = zeros(num_h_modes,num_time_points);

get_h_ddot = @(h,h_dot,W_I,W_C,W_S,w_f) W_I\(w_f - W_C*h_dot - W_S*h);

for iTime = 1:num_time_points
    h_ddot_one(:,iTime) = get_h_ddot(h(:,iTime),h_dot(:,iTime),validation_terms_one{1}(:,:,iTime),validation_terms_one{2}(:,:,iTime),validation_terms_one{3}(:,:,iTime),validation_terms_one{4}(:,iTime));
    h_ddot_two(:,iTime) = get_h_ddot(h(:,iTime),h_dot(:,iTime),validation_terms_two{1}(:,:,iTime),validation_terms_two{2}(:,:,iTime),validation_terms_two{3}(:,:,iTime),validation_terms_two{4}(:,iTime));

    verification_error_all(:,iTime) = abs(h_ddot_one(:,iTime) -  h_ddot_two(:,iTime))./abs( h_ddot_one(:,iTime));

    %%% debug
    % convective_acc_one(:,iTime) = validation_terms_one{1}(:,:,iTime)\validation_terms_one{2}(:,:,iTime)*h_dot(:,iTime);
    % convective_acc_two(:,iTime) = validation_terms_two{1}(:,:,iTime)\validation_terms_two{2}(:,:,iTime)*h_dot(:,iTime);
    % 
    % stiffness_acc_one(:,iTime) = validation_terms_one{1}(:,:,iTime)\validation_terms_one{3}(:,:,iTime)*h(:,iTime);
    % stiffness_acc_two(:,iTime) = validation_terms_two{1}(:,:,iTime)\validation_terms_two{3}(:,:,iTime)*h(:,iTime);
    % 
    % force_acc_one(:,iTime) = validation_terms_one{1}(:,:,iTime)\validation_terms_one{4}(:,iTime);
    % force_acc_two(:,iTime) = validation_terms_two{1}(:,:,iTime)\validation_terms_two{4}(:,iTime);
end


h_ddot_max = max(abs(h_ddot_one),[],2);
largest_h_ddot = max(h_ddot_max);

r_ddot_max = max(abs(r_ddot),[],2);
largest_r_ddot = max(r_ddot_max);

small_acceleration = (h_ddot_max < largest_h_ddot/10 | h_ddot_max < largest_r_ddot/1000);
verification_error_all(small_acceleration,:) = 0;
verification_error = max(max(verification_error_all,[],2));
end