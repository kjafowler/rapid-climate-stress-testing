function [M0, M1, A, B, C] = Matalas_define(Matalas_RawData)

    % first, convert table to array
    RawArray = table2array(Matalas_RawData);

    % lag-0 (M0) population covariance matrix
    lag = 0; M0 = XCovMat (RawArray', lag);
    
    % lag-1 (M1) population covariance matrix
    lag = 1; M1 = XCovMat (RawArray', lag);    
    
    % compute matrix A
    M0I = M0;
    M0I_inv = MatInv_GaussJordan(M0I);
    A = mtimes(M1, M0I_inv); % could also write A = M1*M0I
    
    % compute matrix C
    D1 = mtimes(M1, M0I_inv); % not clear on how this is different to A
    M1T = MatTrans(M1);
    D2 = mtimes(D1,M1T);
    C = M0 - D2; 

    % compute matrix B
    B = GetMatrixB(C); 

end

function OutMatrix = XCovMat (InMatrix, lag)

    % computes covariance matrix between X(t) and X(lag t)
    % I suspect Matlab can do this more natively than the
    % below brute force code... but I want to make it identical
    % where possible to Rory's Fortran.
    
    ns = size(InMatrix, 1); % number of sites
    num = size(InMatrix, 2); % number of data points (years)
    
    for i = 1:ns
        for j = 1:ns
            sij = 0.0;
            ks = 1 + lag;
            for k = ks:num
                sij = sij + InMatrix(i,k) * InMatrix(j,k-lag);
            end
            a(i,j) = sij / num;      
        end
    end
    
    OutMatrix = a; 

end

function OutMatrix = MatInv_GaussJordan(InMatrix)

	% Computes inverse and determinant of matrix A by the Gauss-Jordan method 
    
    a = InMatrix; 
    ns = size(a, 1); ns_check = size(a, 2); 
    if ns~=ns_check, error('Error: input matrix is not square.'); end
    ipiv(1:ns) = 0; 
    
    for i = 1:ns
       big = 0.0;
       for j = 1:ns
          if (ipiv(j)~=1) 
             for k = 1:ns
                if ipiv(k)==0
                    if abs(a(j,k))>=big
                       big = abs(a(j,k));
                       irow = j;
                       icol = k;
                    end
                elseif ipiv(k)>1
                    error('Matrix M0 is singular - inverse cannot be found');
                end
             end
          end
       end
       ipiv(icol) = ipiv(icol) + 1;
       if (irow~=icol)
           for l = 1:ns
              dum = a(irow,l);
              a(irow,l) = a(icol,l);
              a(icol,l) = dum;
           end
       end
       indxr(i) = irow;
       indxc(i) = icol;
       if a(icol,icol)==0.0
           error('Matrix M0 is singular - inverse cannot be found');
       end
       pivinv = 1.0 / a(icol,icol);
       a(icol,icol) = 1.0;
       for l = 1:ns
          a(icol,l) = a(icol,l) * pivinv;
       end
       for ll = 1:ns
          if ll~=icol
             dum = a(ll,icol);
             a(ll,icol) = 0.0;
             for l = 1:ns
                 a(ll,l) = a(ll,l) - a(icol,l) * dum;
             end
          end
       end
    end   
    
    for l = ns:-1:1
       if indxr(l)~=indxc(l)
          for k = 1:ns
             dum = a(k,indxr(l));
             a(k,indxr(l)) = a(k,indxc(l));
             a(k,indxc(l)) = dum;
          end
       end
    end
    
    OutMatrix = a;

end

function OutMatrix = MatTrans(InMatrix) 
    
    % computes transpose matrix (B) of input matrix (A).
    a = InMatrix; 
    
    ns = size(a, 1); ns_check = size(a, 2); 
    if ns~=ns_check, error('Error: input matrix is not square.'); end
    
    for i = 1:ns
       for j = 1:ns
          b(j,i) = a(i,j);
       end
    end
    
    OutMatrix = b; 
end

function OutMatrix = GetMatrixB(InMatrix)
    
    % compute matrix B based on input of matrix C
    % (kf: I legit don't have more info than this)
    
    c = InMatrix;
    ns = size(c, 1); ns_check = size(c, 2); 
    if ns~=ns_check, error('Error: input matrix is not square.'); end    
    
    i = 1;
    j = 1;
    if c(1,1)<=0.0
       error('Unable to determine matrix B');
    end
    b(1,1) = sqrt(c(1,1));
    for i = 2:ns
       b(i,1) = c(i,1) / b(1,1);
       if i~=2
          for j = 2:i-1
             sij = 0.0;
             for k = 1:j-1
                sij = sij + (b(i,k) * b(j,k));
             end
             if b(j,j)<=0.0
                b(i,j) = 0.0;
             else
                b(i,j) = (c(i,j) - sij) / b(j,j);
             end
          end
       end
       sij = 0.0;
       for j = 1:i-1
          sij = sij + (b(i,j) * b(i,j));
       end
       if (c(i,i)-sij)<=0.0
          b(i,i) = 0.0;
       else
          b(i,i) = sqrt(c(i,i) - sij);
       end
    end
    % set zero elements of b
    for i = 1:ns-1
       for j = i+1:ns
          b(i,j) = 0.0;
       end
    end
    
    OutMatrix = b; 
end