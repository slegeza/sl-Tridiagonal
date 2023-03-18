function [x] = tridiagonal(A,B)
    %Our matrix equation is of the form:
    % / b1 c1 0 0 0  \ / x1 \   / r1 \
    % | a2 b2 c2 0 0 | | x2 |   | r2 |
    % | 0 a3 b3 c3 0 | | x3 | = | r3 |
    % | 0 0 a4 b4 c4 | | x4 |   | r4 |
    % \ 0 0 0 a5 b5  / \ x5 /   \ r5 /
    %
    % so I will create arrays for the a, b, c, and r values
    %

    nA = length(A);
    nB = length(B);
    r = B; %to make coding below clearer
    b = diag(A);
    a = diag(A,-1);
    c = diag(A, 1);
    
    %Create empty arrays
    x = [0];
    beta = [0];
    rho = [0];
    
    %Initial values for beta and rho
    beta(1) = b(1);
    rho(1) = r(1);
    
    %Remaining betas and rhos
    for i=2:nB
       %In these equations, we use a(i-1) in place of a(i) in the text,
       %otherwise the index runs out of range; this still works, because
       %the a and c arrays are 1 index shorter than b. For example, b3 and
       %a3 are adjacent on the diagram, and the two values in this program
       %would be b(3) and a(2)

       beta(i) = b(i) - (a(i-1)*c(i-1))/beta(i-1) ;
       rho(i)  = r(i) - (a(i-1)*rho(i-1))/beta(i-1);
    end
    
    %Final x value via eq. 2.22
    x(nB) = rho(nB)/beta(nB);
    
    %Iterate again to solve for x vector
    for i=1:nB-1
       x(nB-i) = (rho(nB-i) - c(nB-i)*x(nB-i+1))/beta(nB-i) ;
    end