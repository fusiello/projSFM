clear all

T = rand(4, 4); 

% transform A with T --> B
for  k = 1:2
     A{k}  = rand(4, 4); 
     B{k}  = A{k} * T  * 10*randn;
end  

T1 = relative_projectivity( A, B );
% find T s.t. B{k} ~ A{k} * T ;

% check that T1 == T
projective_distance({T1},{T},'sphere')


% transform A with T1 --> B1
for  k = 1:2
     B1{k}  = A{k} * T1;
end  
% check that A*T1 == A*T
projective_distance(B1, B)


% check also the call with two return arguments
[T1, B1] = relative_projectivity(A,B);
projective_distance(B1, B)
