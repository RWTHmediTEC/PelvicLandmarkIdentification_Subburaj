function indices = searchAdjM(mRow,rRow )
% Compare spatial relationship (coded in M/L:A/P:I:S sequence) between two
% rows of data.
% mRow contains spatial relationship of desired landmarks, found in R_A.
% rRow contains relationship between candidate and all surface regions.
% Indidices returns a list of all candidates fulfilling spatial
% relationship
%
% AUTHOR: Felix Krooﬂ
% COPYRIGHT (C) 2016 - 2019 Felix Krooﬂ
% LICENSE: EUPL v1.2
%

indices=struct('indices',{});

for j=1:length(mRow)
    x(1) = 0;
    x(2) = 0;
    x(3) = 0;
    split = strsplit(mRow{1,j},':');
    % If spatial relationship is of type 'x', relationship cannot be
    % precisely determined. Therefore candidates relation is checked for
    % both M&L, A&P, or S&I.
    % If there is no 'x' in coded sequence, a regular check for spatial
    % relationship can be made
    for k=1:3
        if split{1,k} =='x'
            x(k) = 1;
        end
    end
        
    if x(1) == 1
        % If x(1) = 'x' M&L has to be checked
       str = ['M' ':' split{1,2} ':' split{1,3}];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = index;
       
       str = ['L' ':' split{1,2} ':' split{1,3}];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = [indices(j).indices index];
       
       str = ['-' ':' split{1,2} ':' split{1,3}];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = [indices(j).indices index];
       
       indices(j).indices = unique(indices(j).indices);
    elseif x(2) == 1
         % If x(2) = 'x' A&P has to be checked
       if x(3) == 1
            % If x(3) = 'x' I&S has to be checked
            str = [split{1,1} ':' 'A' ':' 'S'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = index;
            
            str = [split{1,1} ':' 'A' ':' 'I'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];       
            
            str = [split{1,1} ':' 'A' ':' '-'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];

            str = [split{1,1} ':' 'P' ':' 'S'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' 'P' ':' 'I'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' 'P' ':' '-'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' '-' ':' 'S'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' '-' ':' 'I'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' '-' ':' '-'];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            indices(j).indices = unique(indices(j).indices);
       else
            str = [split{1,1} ':' 'A' ':' split{1,3}];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = index;
            
            str = [split{1,1} ':' 'P' ':' split{1,3}];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            str = [split{1,1} ':' '-' ':' split{1,3}];
            lm = strfind(rRow, str);
            index = find(not(cellfun('isempty', lm)));
            indices(j).indices = [indices(j).indices index];
            
            indices(j).indices = unique(indices(j).indices);
       end
    elseif x(3) == 1
       % If x(3) = 'x' I&S has to be checked
       str = [split{1,1} ':'  split{1,2} ':' 'S'];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = index;
       
       str = [split{1,1} ':'  split{1,2} ':' 'I'];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = [indices(j).indices index];
       
       str = [split{1,1} ':'  split{1,2} ':' '-'];
       lm = strfind(rRow, str);
       index = find(not(cellfun('isempty', lm)));
       indices(j).indices = [indices(j).indices index];
       
       indices(j).indices = unique(indices(j).indices);
    else
        lm = strfind(rRow, mRow{1,j});
        index = find(not(cellfun('isempty', lm)));  
        indices(j).indices = index;
    end
end


end

