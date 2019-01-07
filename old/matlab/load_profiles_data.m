
folder = 'data\\profiles';
extension = '.dat';

directory = dir(folder);

%init struct array
% profiles_data(1).name = '';
% profiles_data(1).x = 0;
% profiles_data(1).y = 0;
% profiles_data(1).m = 0;
% profiles_data(1).p = 0;
% profiles_data(1).camber = 0;
% profiles_data(1).camber_deg = 0;

index = 0;

for name = {directory.name}
    
    % graffe servono perchè è una cell
    if length(name{1})>4 && strcmp(name{1}(end-3: end), extension) && strcmp(name{1}(1:4), 'naca')
        index = index + 1;
        P = load_profile([folder,'\\',name{1}]);
        profiles_data(index)=P;
    end
    
%     if isempty(profiles_data)
%         array_of_struct = repmat(profiles_data(1), length({directory.name})-2, 1 );
%     end
end
