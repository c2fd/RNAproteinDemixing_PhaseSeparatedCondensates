%Code is used to generate phase separation simulation. File will run
%independently from user input or data files. Will output a 2D image 
%of both phase separated species as subplots through time. 
%Code was created by Dr. Jia Zhao for Gasior and Zhao et al. 2017 in PRE. 
%Email: jia.zhao@usu.edu, kgasior@fsu.edu

function  LLSolver()
% Control Parameters
Ny = 128;
Nz = 128;
StableC0 = 10.;

dt    = 1.0e-3;
t_end = 1000;
time  = 0.;
round = 0;

Ly = 1.;
Lz = 1.;

hy = Ly/Ny;
hz = Lz/Nz;

%Model parameters
eps2    = 1.0e-4;
lam_N1  = 5e-1;
lam_R   = 5.8e-1;
lam_W   = 37.5e-1;
lam_N2  = 4.4e-1;
c1 = 1e0;
c2 = 1e-2;
c3 = 1e-1;
c4 = 1e-2;
Ga = 1.5;
alpha = 0.1;
A = 0.0;

%Set initial conditions
N1 = zeros(Ny,Nz);
R = rand(Ny,Nz);
W = 1-R;
N2 = zeros(Ny,Nz);

% set data for old time steps
N1_old = N1;
R_old = R;
W_old = W;
N2_old = N2;




%Use latter for FFT transform
N = Ny; M = Nz;
Leig  = (((2*cos(pi*(0:N-1)'/(N)))-2)*ones(1,M)) +(ones(N,1)*((2*cos(pi*(0:M-1)/(M)))-2));
Leigh= Leig/(hy*hz);
Seig = Leig/(hy*hz);


% 
% start loops
while (time < t_end)
    
    %update time and rounds
    time = time + dt;
    round= round +1;
    
    %extropolation data
    N1_bar = 2*N1 - N1_old;
    R_bar   = 2*R - R_old;
    W_bar   = 2*W - W_old;
    N2_bar   = 2*N2 - N2_old;

%     
    %Solver go-go-go
    N1_np1 = solve_N1_fft(N1,N1_old,N1_bar,N2_bar,W_bar,R_bar);
    W_np1   = solve_W_fft(W,W_old,W_bar,R_bar,N1_bar,N2_bar);
    R_np1   = solve_RNA_fft(R,R_old,R_bar,W_bar,N1_bar);
    N2_np1   = solve_N2_fft(N2,N2_old,N2_bar,N1_bar,W_bar);
    
    %update data
    N1_old = N1; N1 = N1_np1;
    W_old = W; W = W_np1;
    R_old = R; R = R_np1;
    N2_old = N2; N2 = N2_np1;
%     
    if mod(round,1000) == 0
        subplot(1,2,1);image(N1,'CDataMapping','scaled');colorbar('FontSize',12);caxis([0 1]);xlabel('N1','FontSize',15);set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
        subplot(1,2,2);image(N2,'CDataMapping','scaled');colorbar('FontSize',12);caxis([0 1]);xlabel('N2','FontSize',15);set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
%         subplot(2,3,4);image(W,'CDataMapping','scaled');colorbar('FontSize',12);xlabel('P','FontSize',15);set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
%         subplot(2,3,5);image(R,'CDataMapping','scaled');colorbar('FontSize',12);xlabel('R','FontSize',15);set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
        sprintf('Visualize time t=%f',round*dt)
        pause(2);
    end




end


    function N2_np1 = solve_N2_fft(N2,N2_old,N2_bar,N1_bar,W_bar)
        A_N2 = 1.5/dt * ones(N,M) - (lam_N2*StableC0*Leigh) + (lam_N2*eps2*Leigh.*Leigh);
        F_N2 = lam_N2* (chemical_potential_N2(N1_bar,N2_bar)-StableC0*N2_bar);
        hat_rhs = 1./(2*dt)*(4*dct2(N2)-dct2(N2_old)) + Seig.*dct2(F_N2) + dct2(reactive_N2(N2_bar,W_bar,N1_bar));
        hat_N2 = hat_rhs./A_N2;
        N2_np1 = idct2(hat_N2);
    end

    function N1_np1 = solve_N1_fft(N1,N1_old,N1_bar,N2_bar,W_bar,R_bar)
        A_N1 = 1.5/dt * ones(N,M) - (lam_N1*StableC0*Leigh) + (lam_N1*eps2*Leigh.*Leigh);
        F_N1 = lam_N1* (chemical_potential_N1(N1_bar,N2_bar)-StableC0*N1_bar);
        hat_rhs = 1./(2*dt)*(4*dct2(N1)-dct2(N1_old))+ Seig.*dct2(F_N1) + dct2(reactive_N1(N1_bar,W_bar,R_bar,N2_bar));
        hat_N1 = hat_rhs./A_N1;
        N1_np1 = idct2(hat_N1);
    end

    function R_np1 = solve_RNA_fft(R,R_old,R_bar,W_bar,N1_bar)
        %solver for RNA
        A_R = 1.5/dt * ones(N,M) - lam_R * Leigh;
        hat_rhs = 1./(2*dt)*(4*dct2(R)-dct2(R_old)) + dct2(reactive_R(N1_bar,W_bar,R_bar));
        hat_R = hat_rhs./A_R;
        R_np1 = idct2(hat_R);
    end

    function W_np1 = solve_W_fft(W,W_old,W_bar,R_bar,N1_bar,N2_bar)
        %solver for Whi3
        A_W = 1.5/dt * ones(N,M) - lam_W * Leigh;
        hat_rhs = 1./(2*dt)*(4*dct2(W)-dct2(W_old)) + dct2(reactive_W(N1_bar,W_bar,R_bar,N2_bar));
        hat_W = hat_rhs./A_W;
        W_np1 = idct2(hat_W);
    end




function FN1 = chemical_potential_N1(myN1,myN2)
        %chemical potential of N1; delta F/ delta phi
        FN1 = 2.*myN1.*(1.-(myN1 + myN2)).^2 - 2.*(myN1.^2).*(1.-(myN1+myN2))-2.*(myN2.^2).*(1.-(myN1+myN2))...
        +alpha.*( -Ga.*exp(-Ga.*myN1).*(myN1<A)...
        -Ga.*exp(-Ga.*(1+myN1-myN2)).*((1+myN1-myN2)<A)...
        +Ga.*exp(-Ga.*(1+myN2-myN1)).*((1+myN2-myN1)<A));



end

function FN2 = chemical_potential_N2(myN1,myN2)
        %chemical potential of N2 ; delta F/ delta phi
        FN2 = 2.*myN2.*(1.-(myN1 + myN2)).^2 - 2.*(myN1.^2).*(1.-(myN1+myN2))-2.*(myN2.^2).*(1.-(myN1+myN2))...        
        +alpha.*( -Ga.*exp(-Ga.*myN2).*(myN2<A)...
        +Ga.*exp(-Ga.*(1+myN1-myN2)).*((1+myN1-myN2)<A)...
        -Ga.*exp(-Ga.*(1+myN2-myN1)).*((1+myN2-myN1)<A));

end

    function rhs_N2 = reactive_N2(N2_bar,W_bar,N1_bar)
        %reactive terms for N1
        rhs_N2 = c3 *N1_bar.*W_bar - c4 *N2_bar;
    end

    function rhs_N1 = reactive_N1(N1_bar,W_bar,R_bar,N2_bar)
        %reactive terms for N2
        rhs_N1 = c1 * W_bar.* R_bar - c2 * N1_bar-c3*N1_bar.*W_bar +c4*N2_bar;
    end

    function rhs_R = reactive_R(N1_bar,W_bar,R_bar)
        %reactive terms for RNA
        rhs_R = -c1 *W_bar.* R_bar + c2 * N1_bar;
    end

    function rhs_W = reactive_W(N1_bar,W_bar,R_bar,N2_bar)
        %reactive terms for Whi3
        rhs_W = -c1 * W_bar .* R_bar + c2 * N1_bar-c3*N1_bar.*W_bar+c4*N2_bar;
    end



end


