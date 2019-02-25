function [a0, trans_P] = transition(p) % prob of state (i,j) goes to
% state (i', j') dimension: L^2 X L^2
a0_temp = transprob(1,1,p)' * transprob(1,1,p);
a0 = zeros(1, 30^2);
    a0(1,1) = a0_temp(2,2);
    a0(1,2) = a0_temp(2,3);
    a0(1,31) = a0_temp(3,2);
    a0(1,32) = a0_temp(3,3);

trans_P = zeros(30^2, 30^2);
for i = 2:29
    for j = 2:29
        p_temp = transprob(i,j,p)' * transprob(j,i,p);
        trans_P(30 * (i-1) + j, 30 * (i-2) + j-1) = p_temp(1,1);
        trans_P(30 * (i-1) + j, 30 * (i-2) + j) = p_temp(1,2);
        trans_P(30 * (i-1) + j, 30 * (i-2) + j+1) = p_temp(1,3);
        trans_P(30 * (i-1) + j, 30 * (i-1) + j-1) = p_temp(2,1);
        trans_P(30 * (i-1) + j, 30 * (i-1) + j) = p_temp(2,2);
        trans_P(30 * (i-1) + j, 30 * (i-1) + j+1) = p_temp(2,3);
        trans_P(30 * (i-1) + j, 30 * i + j-1) = p_temp(3,1);
        trans_P(30 * (i-1) + j, 30 * i + j) = p_temp(3,2);
        trans_P(30 * (i-1) + j, 30 * i + j+1) = p_temp(3,3);
    end
end

% last state is (1,1)
trans_P(1,1) = a0_temp(2,2);
trans_P(1,2) = a0_temp(2,3);
trans_P(1,31) = a0_temp(3,2);
trans_P(1,32) = a0_temp(3,3);

% last state is (1,2)
a2 = transprob(1,2,p)' * transprob(2,1,p);
trans_P(2,1) = a2(2,1);
trans_P(2,2) = a2(2,2);
trans_P(2,3) = a2(2,3);
trans_P(2,31) = a2(3,1);
trans_P(2,32) = a2(3,2);
trans_P(2,33) = a2(3,3);

% last state is (2,1)
a22 = a2';
trans_P(31,1) = a22(1,2);
trans_P(31,2) = a22(1,3);
trans_P(31,31) = a22(2,2);
trans_P(31,32) = a22(2,3);
trans_P(31,61) = a22(3,2);
trans_P(31,62) = a22(3,3);

% last state is (30,30)
a144 = transprob(30,30,p)' * transprob(30,30,p);
trans_P(30*30, 28*30+29) = a144(1,1);
trans_P(30*30, 28*30+30) = a144(1,2);
trans_P(30*30, 29*30+29) = a144(2,1);
trans_P(30*30, 29*30+30) = a144(2,2);

% last state is (30,29)
a143 = transprob(30,29,p)' * transprob(29,30,p);
trans_P(30*30-1, 28*30+28) = a143(1,1);
trans_P(30*30-1, 28*30+29) = a143(1,2);
trans_P(30*30-1, 28*30+30) = a143(1,3);
trans_P(30*30-1, 29*30+28) = a143(2,1);
trans_P(30*30-1, 29*30+29) = a143(2,2);
trans_P(30*30-1, 29*30+30) = a143(2,3);

% last state is (29,30)
a142 = a143';
trans_P(28*30+30, 27*30+29) = a142(1,1);
trans_P(28*30+30, 27*30+30) = a142(1,2);
trans_P(28*30+30, 28*30+29) = a142(2,1);
trans_P(28*30+30, 28*30+30) = a142(2,2);
trans_P(28*30+30, 29*30+29) = a142(3,1);
trans_P(28*30+30, 29*30+30) = a142(3,2);

end