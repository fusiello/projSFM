function [nullity_est, s_small, stats] = estimate_nullity(J, k, tol, mode)
%ESTIMATE_NULLITY Estimate the numerical nullity of a matrix.
%
%   [nullity_est, s_small, stats] = ESTIMATE_NULLITY(J, k, tol, mode)
%   estimates the numerical nullity of J by computing the k smallest
%   singular values, or approximations of them, and counting how many are
%   below the threshold tol.
%
%   The input matrix J is assumed to have already been row-normalized by
%   the caller. Under this assumption, norm(J,'fro') = sqrt(size(J,1)),
%   which is stored as a Frobenius scale / upper bound for sigma_max. It is
%   not, in general, the actual largest singular value.
%
%   Inputs:
%       J       Matrix whose numerical nullity is to be estimated.
%               Sparse input is supported and expected for eigs-based modes.
%
%       k       Number of smallest singular values/eigenvalues to compute.
%               Default: 1.
%
%       tol     Numerical threshold. Singular values smaller than tol are
%               counted as numerically zero.
%               Default: 5e-8.
%
%       mode    Method used to estimate the smallest singular values:
%
%                   'svd', 'rank', 'small'
%                       Compute the full economy SVD of full(J). This is the
%                       most direct method, but can be expensive and converts
%                       J to dense format.
%
%                   'svds'
%                       Use svds(J,k,'smallest') to compute the k smallest
%                       singular values directly.
%
%                   'eigs_sm', 'large'
%                       Form the normal matrix M = J'*J, symmetrize and
%                       lightly regularize it, then compute the k eigenvalues
%                       of smallest absolute value using eigs.
%
%                   'eigs'
%                       Same normal-matrix approach, but uses shift-invert
%                       eigs with shift 1e-11.
%
%                   'auto'
%                       Try 'smallestabs' first, then fall back to
%                       shift-invert if eigs reports a singular or
%                       ill-conditioned problem.
%
%                   'own'
%                       Custom eigs-based branch using G = J'*J, explicit
%                       ARPACK options, and svds fallback.
%
%               Default: 'eigs_sm'.
%
%   Outputs:
%       nullity_est
%               Estimated numerical nullity:
%
%                   nullity_est = #{ i : s_small(i) < tol }.
%
%               Note that this only counts among the k computed singular
%               values. If k is smaller than the true nullity, the returned
%               value is truncated by k.
%
%       s_small
%               The computed smallest singular values, sorted in ascending
%               order.
%
%       stats   Structure with diagnostic information. Common fields include:
%
%                   stats.kUse
%                   stats.thr
%                   stats.method
%                   stats.nullity_est
%                   stats.s_small
%                   stats.smax
%
%               For eigs-based modes, stats may also contain:
%
%                   stats.eigs_used
%                   stats.delta
%                   stats.ok_eigs
%
%   Remarks:
%       The eigs-based modes operate on J'*J. Therefore singular values are
%       recovered as sqrt(lambda), where lambda is an eigenvalue of J'*J.
%       A small diagonal regularization delta*I is added in those modes to
%       improve numerical behavior. This lifts exact zero singular values to
%       approximately sqrt(delta), so delta should be kept well below tol^2
%       if nullity estimation near tol is important.
%
%       The final value called stats.smax in the current implementation is
%       sqrt(size(J,1)) under the row-normalization assumption. This is a
%       Frobenius-norm scale and an upper bound for sigma_max, not generally
%       the true largest singular value.
%
%   See also SVD, SVDS, EIGS.

if nargin < 2 || isempty(k),   k = 1;       end
if nargin < 3 || isempty(tol), tol = 5e-8; end
if nargin < 4 || isempty(mode), mode = 'eigs_sm'; end

verbose = false;
stats = struct();

[m,p] = size(J);

switch mode

    case {'svd','rank', 'small'}
        disp('Computing the full rank with svd...')
        s = svd(full(J), 'econ');
        s_small = sort(s, 'ascend'); % singular values

        %s_small(1:min(k, numel(s)));
        % smax = s(end);            % largest singular value

        % Matlab's rank automatic threshold
        % th = max(size(J)) * eps(norm(full(J)));
        % nullity_est = sum(s_small < max(th, tol));
        % stats.thr = max(th, tol);

        nullity_est = sum(s_small < tol);

        stats.kUse = numel(s);
        stats.thr = tol;
        stats.method = 'svd';

    case 'svds'
        disp('Computing the smallest singular values with svds...')

        % smax = svds(J,1,'largest','Display',verbose);
        % smax = max(real(smax),0);

        s_small = svds(J,k,'smallest','Display',verbose);
        s_small = sort(max(real(s_small),0),'ascend');

        nullity_est = sum(s_small < tol);

        stats.kUse = k;
        stats.thr = tol;
        stats.method = 'svds';


        case 'gt_nulls'

        disp('Computing the nullspace with gt_nulls...')

        [N,flag, nullity_est] = gt_nulls(J,3);
       
        stats.kUse = Inf;
        stats.thr = tol;
        stats.method = 'gt_nulls';
        stats.method = 'gt_nulls';
        stats.gt_nulls_flag = flag;

        % artificial
        s_small = ones(1,k);
        s_small(1:nullity_est) =  0;


    case {'eigs_sm', 'large', 'eigs', 'auto'}
        disp('Computing the smallest eigenvalue with eigs...')

        % enforce symmetry
        [M, delta] = enforce_symmetry(J);
        assert(issparse(M));

        % largest for normalization
        % lambda_max = eigs(M,1,'largestabs','Display',verbose);
        % lambda_max = max(real(lambda_max),0);
        % smax = sqrt(lambda_max);

        switch mode
            case {'eigs_sm', 'large'}
                wstate = warning;
                warning('off', 'MATLAB:singularMatrix');
                %warning('off', 'MATLAB:eigs:IllConditionedA');

                lambda = eigs(M,k,'smallestabs','Display',verbose);% < -------------
                eigs_used = 'smallestabs';
                warning(wstate);

            case {'eigs'}
                %sigma = delta + 1e-3 * tol^2; % try this some day
                %sigma  = 1e-11;
                sigma = 0; % M is already shifted by delta
                lambda = eigs(M,k,sigma,'Display',verbose);% < -------------
                eigs_used = 'shift-invert';

            case {'auto'}
                % try smallestabs first and fallback to shift-invert

                wstate = warning;

                try
                    warning('error', 'MATLAB:singularMatrix');
                    warning('error', 'MATLAB:eigs:IllConditionedA');

                    try
                        lambda = eigs(M,k,'smallestabs','Display',verbose);
                        eigs_used = 'smallestabs';
                    catch
                        lambda = eigs(M,k,0,'Display',verbose);
                        eigs_used = 'shift-invert';
                    end

                    warning(wstate);

                catch ME
                    warning(wstate);
                    rethrow(ME);
                end
        end

        lambda = sort(max(real(lambda),0),'ascend');
        s_small = sqrt(lambda); % singular values

        nullity_est = sum(s_small < tol);

        stats.kUse = k;
        stats.thr = tol;
        stats.method = 'eigs';
        stats.eigs_used = eigs_used;
        stats.delta = delta;

        % 
        % case 'own'
        %     disp('Computing smallest eigenvalues with custom eigs branch...')
        %     % Original 4dof version
        % 
        %     % ----- large: row-normalize, Gram, eigs(shift=0), fallback svds ------
        % 
        %     % already done in the caller
        %     % (1) Row-normalize J
        %     % rn = sqrt(sum(J.^2,2));
        %     % rn = full(rn);
        %     % rn(rn==0) = 1;
        %     % Dinv = spdiags(1./rn, 0, m, m);
        %     % Jn = Dinv * J;
        % 
        %     % (2) Gram matrix
        %     G = J' * J;
        %     G = 0.5*(G + G'); % enforce symmetry
        % 
        %     trG = full(trace(G));
        %     if ~isfinite(trG) || trG <= 0, trG = 1; end
        %     delta = 1e-22 * (trG / p);
        %     G = G + delta * speye(p);
        % 
        %     % (3) smallest eigs near 0
        %     lam_small = nan(k,1);
        % 
        %     opts = struct();
        %     opts.isreal = isreal(G);
        %     opts.disp   = 0;
        %     opts.p      = max(4*k + 10, 40);
        %     opts.maxit  = 8000;
        %     opts.tol    = 1e-8;
        %     opts.v0     = randn(p,1);
        % 
        %     try
        %         lam_small = eigs(G, k, 0, opts);
        %         ok_eigs = all(isfinite(lam_small));
        %     catch
        %         ok_eigs = false;
        %     end
        % 
        %     % (4) scale estimate / fallback
        %     if ok_eigs
        %         lam_small = sort(real(lam_small), 'ascend');
        %         lam_small(lam_small < 0) = 0;
        %         s_small = sqrt(lam_small);
        % 
        %         % estimate largest eigenvalue -> largest singular value
        %         % try
        %         %     %lam_max = eigs(G, 1, 'largestreal', opts);
        %         %     %lam_max = max(real(lam_max), 0);
        %         %     %s_max_est = sqrt(lam_max);
        %         % catch
        %         %     %s_max_est = normest(Jn);
        %         % end
        % 
        %     else % eigs failed: fallback
        %         try
        %             s_small = svds(J, k, 0);
        %             s_small = sort(full(s_small), 'ascend');
        %         catch
        %             s_small = nan(k,1); %if also svd fails, return nan
        %         end
        %         %s_max_est = normest(Jn);
        %     end
        % 
        % 
        % nullity_est = sum(s_small < tol);
        % 
        % stats.kUse = k;
        % stats.ok_eigs = ok_eigs;
        % stats.thr = tol;
        % stats.method = 'own';
        % stats.delta = delta;

    otherwise
        error('estimate_nullity:UnknownMode', ...
            'Unknown mode "%s".', mode)
end

% Estimate an upper scale from row-normalization.
% If rows are normalized, norm(J,'fro') = sqrt(m).
% This is not necessarily sigma_max, only a scale / upper estimate.
smax = sqrt(m);

fprintf('Smallest sigma = %.1e\n',s_small(1))
fprintf('Largest  sigma ~ %.1e\n',smax)

stats.nullity_est = nullity_est;
stats.s_small = s_small;
stats.smax    = smax;

end

function [M, delta] = enforce_symmetry(J)
%ENFORCE_SYMMETRY Form a symmetric regularized normal matrix M = J'J.

[m, p] = size(J); %#ok<ASGLU>

M = J' * J;
M = 0.5 * (M + M');

trG = full(trace(M));
if ~isfinite(trG) || trG <= 0
    trG = 1;
end

delta = 1e-22 * (trG / p);
M = M + delta * speye(p);

end