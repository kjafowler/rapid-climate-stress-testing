function Matalas_RawOutput = Matalas_apply(A, B, TotalYrs, info)
    
    % Generate a single replicate using the Matalas method, based on
    % pre-calculated matrices.  
    
    LenRep = TotalYrs;
    
    NumSite = size(A, 1); ns_check = size(A, 2); 
    if NumSite~=ns_check, error('Error: input matrix is not square.'); end    
    
    X1(1:NumSite) = 0.0;
    X2(1:NumSite) = 0.0;
    for i = 1:LenRep
       rn = rand(1, 8); 
       for j = 1:NumSite
          xrn = rn(j);
          % Zstd(j,i) = ppnd(xrn); % tried and failed to convert function ppnd from Fortran to Matlab (see below)
          Zstd(j,i) = norminv(xrn);
          sij = 0.0;
          for k = 1:NumSite
             sij = sij + A(j,k) * X1(k);
          end
          for k = 1:j
             sij = sij + B(j,k)*Zstd(k,i);
          end
          X2(j) = sij;     
       end
       
        tf = isreal(X2);
        if ~tf
            test01 = 1; 
        end       
       
       for j = 1:NumSite
          X1(j) = X2(j);
          V(j,i) = X2(j);
          % if (flag.eq.'N') then
          %    V(j,i) = X2(j);
          % elseif (flag.eq.'S') then
          %    V(j,i) = sdt(j) * X2(j)
          % elseif (flag.eq.'L') then
          %    V(j,i) = av(j) + sdt(j) * X2(j)
          % endif
       end
    end
    
    tf = isreal(V);
    if ~tf
        test01 = 1; 
    end    
    
    % table format for output
    Matalas_RawOutput = array2table(V', 'VariableNames', info.SubareaList);
    
end

% function NormalDeviate = ppnd(P)
%     
%     % Created 04/10/2019 by Keirnan Fowler, University of Melbourne
%     % I'm sure there's an equivalent procedure in Matlab but I want this
%     % to be identical to Fortran code provided by Rory Nathan.
%     
%     % ORIGINAL CODE
%     
% %       REAL FUNCTION PPND(P,IFAULT)
% % !
% % !     PRODUCES NORMAL DEVIATE CORRESPONDING TO LOWER TAIL AREA OF P.
% % !     ALGORITHM AS 111 APPLIED STATISTICS (1977) VOL. 26, P 118
% % !
% %       REAL ZERO, SPLIT, HALF, ONE, A0, A1, A2, A3, B1, B2, B3, B4
% %       REAL C0, C1, C2, C3, D1, D2, P, Q, R, ZABS, ZLOG, ZSQRT
% % !
% %       DATA ZERO, HALF, ONE, SPLIT /0.0E0, 0.5E0, 1.0E0, 0.42E0/
% %       DATA A0 /        2.50662 82388 4E0/, &
% %            A1 /      -18.61500 06252 9E0/, &
% %            A2 /       41.39119 77353 4E0/, &
% %            A3 /      -25.44106 04963 7E0/, &
% %            B1 /       -8.47351 09309 0E0/, &
% %            B2 /       23.08336 74374 3E0/, &
% %            B3 /      -21.06224 10182 6E0/, &
% %            B4 /        3.13082 90983 3E0/
% % !
% % !        HASH SUM AB 143.70383 55807 6
% % !
% %       DATA C0 /       -2.78718 93113 8E0/, &
% %            C1 /       -2.29796 47913 4E0/, &
% %            C2 /        4.85014 12713 5E0/, &
% %            C3 /        2.32121 27685 8E0/, &
% %            D1 /        3.54388 92476 2E0/, &
% %            D2 /        1.63706 78189 7E0/
% % !
% % !       HASH SUM CD 17.43746 52092 4
% % !
% %       ZABS(P) = ABS(P)
% %       ZLOG(P) = ALOG(P)
% %       ZSQRT(P) = SQRT(P)
% % !
% %       IFAULT = 0
% %       Q = P - HALF
% %       IF (ZABS(Q) .GT. SPLIT) GOTO 1
% %       R = Q * Q
% %       PPND = Q * (((A3 * R + A2) * R + A1) * R + A0) / &
% %         ((((B4 * R + B3) * R + B2) * R + B1) * R + ONE)
% %       RETURN
% %    1  R = P
% %       IF (Q .GT. ZERO) R = ONE - P
% %       IF (R .LT. ZERO) GOTO 2
% %       R = ZSQRT(-ZLOG(R))
% %       PPND = (((C3 * R + C2) * R + C1) * R + C0) / ((D2 * R + D1) * R + ONE)
% %       IF (Q .LT. ZERO) PPND = -PPND
% %       RETURN
% %    2  IFAULT = 1
% %       PPND = ZERO
% %       RETURN
% %       END    
%     
%     split = 0.42; 
%     A0 = [  2.50662 82388 4E0];
%     A1 = [-18.61500 06252 9E0];
%     A2 = [ 41.39119 77353 4E0];
%     A3 = [-25.44106 04963 7E0];
%     B1 = [ -8.47351 09309 0E0];
%     B2 = [ 23.08336 74374 3E0];
%     B3 = [-21.06224 10182 6E0];
%     B4 = [  3.13082 90983 3E0];
%     
%     % HASH SUM AB 143.70383 55807 6
%     
%     C0 = [-2.78718 93113 8E0];
%     C1 = [-2.29796 47913 4E0];
%     C2 = [ 4.85014 12713 5E0];
%     C3 = [ 2.32121 27685 8E0];
%     D1 = [ 3.54388 92476 2E0];
%     D2 = [ 1.63706 78189 7E0];
%     
%     % HASH SUM CD 17.43746 52092 4
%     
%     % ZABS(P) = abs(P);
%     % ZLOG(P) = log(P); 
%     % ZSQRT(P) = sqrt(P);
%     
%     Q = P - 0.5;
%     if abs(Q) > split
%         R = P;
%         if Q > 0.0, R = 1 - P; end
%         if R < 0.0 
%             error('IFAULT = 1'); 
%         else
%             R = sqrt(-log(R)); 
%             PPND = (((C3 * R + C2) * R + C1) * R + C0) / ((D2 * R + D1) * R + 1.0); 
%             if Q < 0, PPND = -PPND; end
%         end
%     else
%         R = Q * Q;
%         PPND = Q * (((A3 * R + A2) * R + A1) * R + A0) / ((((B4 * R + B3) * R + B2) * R + B1) * R + 1.0);
%     end
%     
%     NormalDeviate = PPND; 
%     
% end
