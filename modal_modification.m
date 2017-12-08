%% Clear previous iterations and save default workspace from HW 1.

% This script was designed for cell mode, but this block must be run once
% to generate the environmental variables that the rest of the cells use
% for their calculations.
clc; clear; close all

wrDef = [40.6, 79.6, 100, 150.1, 298, 451, 612].*2*pi;
zed = [0.0786 0.0345 0.0516 0.0507 0.0415 0.0269 0.0272];
med = [1 1 1 1 1 1 1]./7;
zetaDef = diag(zed);
mDef = diag(med);
cDef = [2*zetaDef(1,1)*wrDef(1)*mDef(1,1) 0 0 0 0 0 0; 0 2*zetaDef(2,2)*wrDef(2)*mDef(2,2) 0 0 0 0 0; 0 0 2*zetaDef(3,3)*wrDef(3)*mDef(3,3) 0 0 0 0; 0 0 0 2*zetaDef(4,4)*wrDef(4)*mDef(4,4) 0 0 0; 0 0 0 0 2*zetaDef(5,5)*wrDef(5)*mDef(5,5) 0 0; 0 0 0 0 0 2*zetaDef(6,6)*wrDef(6)*mDef(6,6) 0; 0 0 0 0 0 0 2*zetaDef(7,7)*wrDef(7)*mDef(7,7)];
cDef = (diag(cDef)).';
uDef = ones(7);
uTDef = uDef.';
ked = (1/7).*[wrDef(1)^2 wrDef(2)^2 wrDef(3)^2 wrDef(4)^2 wrDef(5)^2 wrDef(6)^2 wrDef(7)^2];
kDef = diag(ked);
save('hw3');

%% Problem 0: Unmodified

% Load environmental variables for cell-mode users.
load('hw3');

% Problem 0 is the same as the cell above; this cell is to ensure all parts
% have been completed.

%% Problem 1: Mass Modification.
% 1-kg weights added to each mass in the system, no change to stiffness.

% Load environmental variables for cell-mode users.
load('hw3');

% Define and quantify changes to the physical system.
delM = [0 0 0 7.93785 7.93785+10.20585 10.20585 0];
dm1 = diag(delM);
Mm1 = mDef + uTDef*dm1*uDef;
Km1 = kDef;

% Run eigenvalue solver and sort output.
[u1Pr,dPr1]=eig(Km1, -Mm1);
% u1Pr = u1Pr(:, [1 3 2]);
% dPr1 = dPr1(:, [1 3 2]);

% Extract modified natural frequencies from eig() output.
modWnRad1 = sqrt(abs(dPr1));
% modWnRad1 = sort(diag(modWnRad1));
modWnHz1 = (modWnRad1/(2*pi)).';

% Scale the [U'] matrix by unity in the first row.
for r = 1:3
    u1PrScaled(:,r) = u1Pr(:,r)/u1Pr(1,r);
end

% Calculate modal mass and modified stiffness matrix.

modMass1 = u1PrScaled.'*Mm1*u1PrScaled;
modStiff1 = u1PrScaled.'*Km1*u1PrScaled;

% Calculate U1: [T] = [U][U']
u1US = uDef*u1PrScaled;
for r = 1:3
    U1(:,r) = u1US(:,r)/u1US(1,r);
end
 
% for z = 1:3
%     zeta1(z)=(cDef(z)./(2.*modWnRad1(z).*modMass1(z,z)));
% end

for z = 1:3
    zeta1(z)=(cDef(z)./(2.*sqrt(modMass1(z,z).*modStiff1(z,z))));
end


%% Problem 2: Stiffness Modification
% A 40MN/m spring is connected between masses m2 and m3.

% Load environmental variables for cell-mode users.
load('hw3');

% Define and quantify changes to the physical system.
dk2 = [0 0 0; 0 40e6 -40e6; 0 -40e6 40e6];
Mm2 = mDef;
Km2 = kDef + uTDef*dk2*uDef;

% Run eigenvalue solver.
[u2Pr,dPr2]=eig(Km2, -Mm2);

% Extract modified natural frequencies from eig() output.
modWnRad2 = sqrt(abs(dPr2));
modWnRad2 = sort(diag(modWnRad2));
modWnHz2 = (modWnRad2/(2*pi))';

% Scale the matrix by the first non-zero element.
u2PrScaled = abs(u2Pr);

% Calculate modal mass and modified stiffness matrix.

modMass2 = u2PrScaled.'*Mm2*u2PrScaled;
modStiff2 = u2PrScaled.'*Km2*u2PrScaled;

% Calculate U1: [T] = [U][U']
u2US = uDef*u2PrScaled;
for r = 1:3
    U2(:,r) = u2US(:,r)/u2US(1,r);
end

for z = 1:3
    zeta2(z)=(cDef(z)./(2.*sqrt(modMass2(z,z).*modStiff2(z,z))));
end
%% Problem 3: Tuned Vibration Absorber
% A 0.5kg mass is attached to m2 by a 1.78MN/m spring.

% Load environmental variables for cell-mode users.
load('hw3');

% Define and quantify changes to the physical system.
dm3 = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0.5];
dk3 = [0 0 0 0; 0 1.78e6 0 -1.78e6; 0 0 0 0; 0 -1.78e6 0 1.78e6];
mDef3 = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 0];
kDef3 = [wrDef(1)^2 0 0 0; 0 2*wrDef(2)^2 0 0; 0 0 3*wrDef(3)^2 0; 0 0 0 0];
uDef3 = [1 1 1 0; 1 -1 -1 0; 1 -1 1 0; 0 0 0 1];
uTDef3 = uDef3.';
Mm3 = mDef3 + uTDef3*dm3*uDef3;
Km3 = kDef3 + uTDef3*dk3*uDef3;
 
cDef3 = [2*zetaDef(1,1)*wrDef(1)*mDef(1,1) 0 0 0; 0 2*zetaDef(2,2)*wrDef(2)*mDef(2,2) 0 0; 0 0 2*zetaDef(3,3)*wrDef(3)*mDef(3,3) 0; 0 0 0 0];
cDef3 = (diag(cDef3)).';        % Check to see how damping ratio is determined for m_a

% Run eigenvalue solver and sort output.
[u3Pr,dPr3]=eig(Km3, -Mm3);
u3Pr = u3Pr(:, [3 4 2 1]);
dPr3 = dPr3(:, [3 4 2 1 ]);

% Extract modified natural frequencies from eig() output.
modWnRad3 = sqrt(abs(dPr3));
modWnHz3 = (modWnRad3/(2*pi))';

% Scale the [U'] matrix by unity in the first row.
for r = 1:4
    u3PrScaled(:,r) = u3Pr(:,r)/u3Pr(1,r);
end

% Calculate modal mass and modified stiffness matrix.

modMass3 = u3PrScaled.'*Mm3*u3PrScaled;
modStiff3 = u3PrScaled.'*Km3*u3PrScaled;

% Calculate U1: [T] = [U][U']
u3US = uDef3*u3PrScaled;
for r = 1:4
    U3(:,r) = u3US(:,r)/u3US(1,r);
end

for z = 1:4
    zeta3(z)=(cDef3(z)./(2.*sqrt(modMass3(z,z).*modStiff3(z,z))));
end
%% Wrapping Up: Data Aggregation and Plotting

% The following line cleans up by deleting the data file (HW1.mat) from 
% your hard drive without confirmation. Comment or remove the line if this
% is not desired.
delete hw3.mat

% Plot Mode 1 shape.
x1 = [1:3];
x2 = [1:4];
figure(1)
plot(x1,uDef(:,1),x1,U1(:,1),x1,U2(:,1),x2,U3(:,1))
title('Mode 1')
xlabel('Degree of Freedom')
ylabel('Relative Motion')
legend('0. Unmodified', '1. Mass Modification', '2. Stiffness Modification', '3. TVA')
axis([1 4 0.5 1.5]);
grid on;


% Plot Mode 2 shape.
x1 = [1:3];
x2 = [1:4];
figure(2)
plot(x1,uDef(:,2),x1,U1(:,2),x1,U2(:,2),x2,U3(:,2))
title('Mode 2')
xlabel('Degree of Freedom')
ylabel('Relative Motion')
legend('0. Unmodified', '1. Mass Modification', '2. Stiffness Modification', '3. TVA')
axis([1 4 -1.5 1.5]);
grid on;


% Plot Mode 3 shape.
x1=[1:3];
x2 = [1:4];
figure(3)
plot(x1,uDef(:,3),x1,U1(:,3),x1,U2(:,3),x2,U3(:,3))
title('Mode 3')
xlabel('Degree of Freedom')
ylabel('Relative Motion')
legend('0. Unmodified', '1. Mass Modification', '2. Stiffness Modification', '3. TVA')
axis([1 4 -3.5 3.5]);
grid on;

% Plot the odd TVA mode shape.
x = [1:4];
y = (U3(:,4)).';
figure(4)
plot(x,y)
title('TVA Mode 4')
xlabel('Degree of Freedom')
ylabel('Relative Motion')
grid on
