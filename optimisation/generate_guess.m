function [c0, param, phases] = generate_guess (param, phases)

%% initialize normalized bounds
UB = ones(1,numel(param.UB));
LB = 0*UB;

%% define population size (default = 1)
popsz = 1;
for i=1:size(param.fg,1)
    if strcmp(param.fg(i,1),'n_individuals')
        popsz=cell2mat(param.fg(i,2));
    end
end

%% load population seed (first individual, default = mean between bounds)
c0 = 0.5*UB;
for i=1:size(param.fg,1)
    if strcmp(param.fg(i,1),'seed')
        if strcmp(param.fg(i,2),'latest')
            cd('./results');
            d = dir('*.mat');
            [dx,dx] = sort([d.datenum]);
            newest = d(dx(end)).name;
            load(newest,'c0')
            cd('../');
        elseif iscellstr(param.fg(i,2))
            load(param.fg{i,2},'c0')
        elseif isnumeric(param.fg(i,2))
            if length(c0) == length(param.fg(i,2))
                c0 = cell2mat(param.fg(i,2));
            else
                disp(['Size of the array loaded: ',num2str(length(param.fg(i,2))),' elements'])
                disp(['Size expected: ',num2str(length(c0)),' elements'])
                error('Verify the argument assigned to the seed in param.fg')
            end
        end
    end
end

%% define rest of the population distribution (default = single individual)
for i=1:size(param.fg,1)
    if strcmp(param.fg(i,1),'distribution')
        if strcmp(param.fg(i,2),'lhs')
            c0 = lhsdesign(popsz,length(UB));
        elseif strcmp(param.fg(i,2),'random')
            c0 = rand(popsz,length(UB));
        end
    end
end

%% apply noise (default = no noise)
for i=1:size(param.fg,1)
    if strcmp(param.fg(i,1),'noise')
        if size(c0,1)==1
            noise = cell2mat(param.fg(i,2));
            c0 = repmat(c0,popsz,1)-0.5*noise+noise*rand(popsz,length(UB));
        else
            disp('Noise can be applied only to a provided guess, not to a distribution')
            disp('either disable noise or distribution fileds in param.fg')
        end
    end
end

%% final check
if size(c0,1)~=popsz
    error ('Check param.fg settings, the population size is not compatible')
end