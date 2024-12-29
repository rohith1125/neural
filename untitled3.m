% Define the patterns P, Q, and B
P = [1 -1 1 -1];
Q = [-1 -1 1 -1];
B = [0 -1 1 0];

% Part a: Create a Hopfield network H that stores the patterns P and Q
H = memstor([P; Q]);
fprintf(H)
% Part b: Calculate the goodness values for H with each of the following patterns: P, Q, and B
goodness_P = calc_specific_goodness(H, P);
goodness_Q = calc_specific_goodness(H, Q);
goodness_B = calc_specific_goodness(H, B);

% Part c: Create a Hopfield network K that stores the patterns P, Q, and B
K = memstor([P; Q; B]);

% Calculate the goodness values for all 3 patterns on network K
goodness_P_K = calc_specific_goodness(K, P);
goodness_Q_K = calc_specific_goodness(K, Q);
goodness_B_K = calc_specific_goodness(K, B);

% Display the results
fprintf('Goodness values for network H:\n');
fprintf('Pattern P: %f\n', goodness_P);
fprintf('Pattern Q: %f\n', goodness_Q);
fprintf('Pattern B: %f\n\n', goodness_B);

fprintf('Goodness values for network K:\n');
fprintf('Pattern P: %f\n', goodness_P_K);
fprintf('Pattern Q: %f\n', goodness_Q_K);
fprintf('Pattern B: %f\n', goodness_B_K);

% Your provided functions:
function mem = memstor(pats)
    % each row of the matrix pats is a pattern
    [np, nd] = size(pats);
    mem = zeros(nd);
    for i = 1:nd
        for j = 1:nd
            if (i ~= j)
                for k = 1:np
                    mem(i,j) = mem(i,j) + pats(k,i) * pats(k,j);
                end
            end
        end
    end
end

function gvals = goodness(hopnet)
    % calculates goodness for all patterns in a Hopfield Network
    gvals = [];
    pmat = [];
    netsize = size(hopnet,1);
    for k = 0:(2^netsize-1)
        pvec = 2*de2bi(k,netsize,'left-msb')-1;
        pmat = [pmat; pvec];
        g = 0;
        for i = 1:(netsize-1)
            for j = (i+1):netsize
                g = g + hopnet(i,j) * pvec(i) * pvec(j);
            end
        end
        gvals = [gvals, g];
    end
    gvals = gvals';
    [pmat, gvals]
end

% New function based on your requirements:
function gval = calc_specific_goodness(hopnet, pattern)
    % Calculate the goodness for a specific pattern
    netsize = length(pattern);
    gval = 0;
    for i = 1:(netsize-1)
        for j = (i+1):netsize
            gval = gval + hopnet(i,j) * pattern(i) * pattern(j);
        end
    end
end

