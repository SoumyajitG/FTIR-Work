function [x] = OMP( A, b, k )

LARGESCALE  = false;
Af  = @(x) A*x;
At  = @(x) A'*x;


% -- Intitialize --
% start at x = 0, so r = b - A*x = b
r           = b;
Ar          = At(r);
N           = size(Ar,1);       % number of atoms
M           = size(r,1);        % size of atoms
unitVector  = zeros(N,1);
x           = zeros(N,1);

indx_set    = zeros(k,1);
indx_set_sorted     = zeros(k,1);
A_T         = zeros(M,k);
A_T_nonorth = zeros(M,k);
slowMode = 0;


for kk = 1:k
    
    % -- Step 1: find new index and atom to add
    [dummy,ind_new]     = max(abs(Ar));
    % Check if this index is already in
%     if ismember( ind_new, indx_set_sorted(1:kk-1) )
%         disp('Shouldn''t happen... entering debug');
%         keyboard
%     end
    
    
    indx_set(kk)    = ind_new;
    indx_set_sorted(1:kk)   = sort( indx_set(1:kk) );
    
    if LARGESCALE
        unitVector(ind_new)     = 1;
        atom_new                = Af( unitVector );
        unitVector(ind_new)     = 0;
    else
        atom_new    = A(:,ind_new);
    end
    
    A_T_nonorth(:,kk)   = atom_new;     % before orthogonalizing and such
    
    
    
    % -- Step 2: update residual
    
    if slowMode
        % The straightforward way:
        x_T = A_T_nonorth(:,1:kk)\b;
        
        % or, use QR decomposition:
%         if kk < 10
% %             [Q,R] = qr( A_T_nonorth(:,1:kk), 0 );
%             [Q,R] = qr( A_T_nonorth(:,1:kk)); % need full "Q" matrix to use "qrinsert"
%             %  For this reason, "qrinsert" is not efficient
%         else
%             % from now on, we use the old QR to update the new one
%             [Q,R] = qrinsert( Q, R, kk, atom_new );
%         end
%         x_T = R\(R'\(A_T_nonorth(:,1:kk)'*b));
        
        
        x( indx_set(1:kk) )   = x_T;
        r   = b - A_T_nonorth(:,1:kk)*x_T;
    else
    
        % First, orthogonalize 'atom_new' against all previous atoms
        % We use MGS
        for j = 1:(kk-1)
%             atom_new    = atom_new - (atom_new'*A_T(:,j))*A_T(:,j);
            % Thanks to Noam Wagner for spotting this bug. The above line
            % is wrong when the data is complex. Use this:
            atom_new    = atom_new - (A_T(:,j)'*atom_new)*A_T(:,j);
        end
        % Second, normalize:
        atom_new        = atom_new/norm(atom_new);
        A_T(:,kk)       = atom_new;
        % Third, solve least-squares problem (which is now very easy
        %   since A_T(:,1:kk) is orthogonal )
        x_T     = A_T(:,1:kk)'*b;
        x( indx_set(1:kk) )   = x_T;      % note: indx_set is guaranteed to never shrink
        % Fourth, update residual:
        %     r       = b - Af(x); % wrong!
        r       = b - A_T(:,1:kk)*x_T;
        
        % N.B. This err is unreliable, since this "x" is not the same
        %   (since it relies on A_T, which is the orthogonalized version).
    end
    
 
    
    if kk < k
        Ar  = At(r); % prepare for next round
    end
    
end
if ~slowMode  % (in slowMode, we already have this info)
 % For the last iteration, we need to do this without orthogonalizing A
 % so that the x coefficients match what is expected.
 x_T =pinv( A_T_nonorth(:,1:kk))*b;
 x( indx_set(1:kk) )   = x_T;
end
end