function TETA      = TETA_BIG(q_)
% q_=[q1 q2 q3 q4]'=[q q4]';
TETA=[q_(4)*eye(3)+a2A(q_(1:3));-q_(1:3)'];