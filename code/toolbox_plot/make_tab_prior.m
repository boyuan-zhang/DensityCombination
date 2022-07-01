function  tab1 = make_tab_prior(prior, paraname)
% make table for prior in cell format

% size info
npara = size(prior,1);

% if paraname is not in the function, them make it empty
if nargin == 1
    paraname = cell(npara,1);
end

% labeling
shape_set = {'Beta', 'Gamma', 'Normal', 'InvGamma', 'Uniform', 'No prior', 'Fixed'};

% make table
prior_tab_head = {'Parameter', 'Distribution', 'Para(1)', 'Para(2)'};
prior_tab = num2cell(prior);


for i=1:1:npara

    % check fixed parameter
    if prior(i,4) == 1
        prior(i,1) = 7; % set prior dist. to "fixed"
        prior_tab{i,2} = prior(i,5); % set Para(1) to fixed value
        prior_tab{i,3} = 'N/A'; % set Para(1) to fixed value
    end
    
    % name of distribution
    pdf_sh = prior(i,1);
    prior_tab{i,1} = shape_set{pdf_sh};
    
end

% make table
tab1 = [prior_tab_head; [paraname prior_tab(:,1:3)]];
