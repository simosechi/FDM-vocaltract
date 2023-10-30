function [] = main(testname)

close all

%%%%%%%%%% Prendere i dati da F_Dati %%%%%%%%%%%%

addpath("Utils\");
dati = F_dati(testname);

%%%%%%%%%%% Costruzione della mesh %%%%%%%%%%%%%%

x_dom = linspace(0,dati.length, dati.n); %dominio nello spazio
t_dom = (0:dati.dt:dati.T-dati.dt); %dominio nel tempo
x = linspace(0,dati.length,2*dati.n).';
t = (0:dati.dt/2:dati.T-dati.dt/2);

%%%%%%%%%%%% VOWELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_x = [0.0025 0.0023 0.002 0.006 0.0105 0.009 0.008 0.0065 0.005 0.0035 0.0027 0.002 0.0015 0.001 0.0015 0.002 0.0017 0.0015 0.0017 0.002 0.0027 0.0035 0.0063 0.008 0.01 0.0135 0.0135 0.0135 0.012 0.011 0.0085 0.006 0.004 0.002 0.0012 0.0005 0.001 0.001];
A = zeros(size(x,1), size(t,2));
vowels = zeros(5,38);
vowels(1:end-1,1:2:end) = [
    0.0025	0.001	0.0038	0.002	0.0008	0.0008	0.0017	0.002	0.0025	0.0018	0.0052	0.008	0.008	0.008	0.008	0.007	0.005	0.005	0.005;     %a
    0.0025  0.0015  0.0065  0.0065  0.008   0.009   0.0105  0.0105  0.0065  0.008   0.005 0.004 0.0032 0.0028 0.0028 0.005 0.008 0.008 0.008;  %e
    0.003 0.002 0.0078 0.0105 0.0105 0.0105 0.0105 0.0105 0.008 0.0042 0.0016 0.0005 0.0005 0.0005 0.001 0.0019 0.004 0.004 0.004;  %i
    0.0025 0.001 0.0048 0.004 0.0022 0.0015 0.001 0.001 0.0022 0.0018 0.004 0.005 0.0068 0.008 0.0105 0.0158 0.013 0.006 0.0037; %o
    ];
vowels(1:end-1,2:2:end-1) = (vowels(1:end-1, 1:2:end-3) + vowels(1:end-1, 3:2:end-1))/2;
vowels(1:end-1, end) = vowels(1:end-1, end-1);
vowels(end, :)= A_x;

couples = minDispVowelsPattern(vowels);

vows = ['A', 'E', 'I', 'O', 'U'];
vowIdx = [1, 2, 3, 4, 5];
pattern = perms(vowIdx);
patternCost = zeros(size(pattern,1),1);

for i = 1:size(pattern,1)
    patternCost(i,1) = couples(pattern(i,1), pattern(i,2)) + couples(pattern(i,2), pattern(i,3)) + couples(pattern(i,3), pattern(i,4)) + couples(pattern(i,4), pattern(i,5));
end


minCost = min(patternCost(patternCost>0), [], "all");
minPattern = find(patternCost == minCost);
for i = 1:size(minPattern,1)
    i1 = pattern(minPattern(i,1),1);
    i2 = pattern(minPattern(i,1),2);
    i3 = pattern(minPattern(i,1),3);
    i4 = pattern(minPattern(i,1),4);
    i5 = pattern(minPattern(i,1),5);
    fprintf(['\n','min pattern ', i, '= ',vows(i1), vows(i2), vows(i3), vows(i4), vows(i5) , '\n']);
end

figure('Name',['Vocal Tract area function for A-E-I vowels']);
for vow=1:length(vows)-2
    stairs(x,vowels(vow, :));
    hold on;
end
legend (['vowel ', vows(1)],['vowel ', vows(2)],['vowel ', vows(3)]);
xlabel('Distance from the glottis [$m$]', 'Interpreter','latex');
ylabel('Area of the cross-section [$m^2$]', 'Interpreter','latex');
hold off;

figure('Name',['Vocal Tract area function for O-U vowels']);
for vow=4:5
    stairs(x,vowels(vow, :));
    hold on;
end
legend (['vowel ', vows(4)],['vowel ', vows(5)]);
xlabel('Distance from the glottis [$m$]', 'Interpreter','latex');
ylabel('Area of the cross-section [$m^2$]', 'Interpreter','latex');
hold off;



Xsize = size(x_dom, 2);
Tsize = size(t_dom, 2);

% impostare valori iniziali



% impulso iniziale dato alla velocity nel tempo  


Gp = glottalpulse(dati.dt, dati.t0, dati.t1, dati.t2, dati.Ag);
Gpsize = size(Gp,2);




CFL = dati.c * dati.dt * dati.n / dati.length;
fprintf(['CFL number   ', num2str(CFL)], '\n');
Vout = zeros(size(vowels,1), Tsize);
L = Tsize;
f = dati.fs*(0:(L/2))/L;
H = zeros(size(vowels,1), length(f));
H_log = zeros(size(vowels,1), length(f));
formants = zeros(5,3);



%%%%%%%%%%%%%%% START LOOP %%%%%%%%%%%%%%%%%%%

for vow = 1:size(vowels, 1)

    for n=1:size(t, 2)
        A(:,n)= vowels(vow, :);
    end


    u = zeros(Xsize, Tsize);
%     u(1,1) = Gp(1,1);
    u(1,1) = 0;
    u(1,3) = 0.08;
    Asize_x = size(A,1);
    Asize_t = size(A,2);    
    p = zeros(Xsize, Tsize);
    Pdelta = zeros(Xsize - 1,1);   
    u_old = u(2:Xsize,1);
    p_old = p(1:Xsize-1,1);
    A_tilde_old = (1./A(2:2:Asize_x-2,2) + 2./A(3:2:Asize_x-1,2) + 1./A(4:2:Asize_x,2))./4;
    A_cap_old = (A(1:2:Asize_x-3,1) + 2.*A(2:2:Asize_x-2,1) + A(3:2:Asize_x-1, 1))./4;
    w_old = 0;

    for i = 2:Tsize
    
%         j = mod(i-1,Gpsize)+1;
%         u(1,i) = Gp(1,j);
    
        L_R = 8/(3*pi.*sqrt(pi.*A(end, 2*i)));
        R_R = 128*dati.c/(9*A(end,2*i)*pi^2);    
        w_new = w_old + dati.dt.*dati.c./L_R.*p_old(Xsize-1);
        u_last = dati.c*p_old(Xsize-1)/R_R + w_new;
        %u_last = 0;
        w_old = w_new;
    
        A_tilde_new = (1./A(2:2:Asize_x-2,i*2) + 2./A(3:2:Asize_x-1,i*2) + 1./A(4:2:Asize_x,i*2))./4;
        D_cap = 0.002.*(A(2:2:Asize_x-2, i*2).^(-3/2) + 2.*A(3:2:Asize_x-1,i*2).^(-3/2)+A(4:2:Asize_x,i*2).^(-3/2))./4;
        d_cap = 1.6.*(A(2:2:Asize_x-2, i*2).^(-3/2) + 2.*A(3:2:Asize_x-1,i*2).^(-3/2)+A(4:2:Asize_x,i*2).^(-3/2))./4;
        a = [1 ; A_tilde_new.*(dati.length./dati.n).^2 + d_cap.*dati.dt.*(dati.length./dati.n).^2 + 2.*D_cap.*dati.dt];
        b = [0; -D_cap(1:Xsize-2).*dati.dt];
        c = -D_cap(1:Xsize-1).*dati.dt;
        M = diag(a) + diag(b,1) + diag(c,-1);
        M(end, end-1:end) = [0,1];
        y = [u(1,i) ; 
            -dati.c.*Pdelta(1:end-1).*dati.length./dati.n.*dati.dt + u_old(1:end-1).*A_tilde_old(1:end-1).*(dati.length./dati.n).^2;
            u_last
            ];
        u_new = linsolve(M,y);  % risolve la matrice M u_new = y ,   u_new = M^-1 * y
    
        
    
        u(2:Xsize,i) = u_new(2:Xsize);
        A_tilde_old = A_tilde_new;    
        Udelta = u(2:Xsize,i) - u(1:Xsize-1,i);
        u_old = u_new(2:Xsize);
    
        A_cap_new = (A(1:2:Asize_x-3,2*i-1) + 2.*A(2:2:Asize_x-2,2*i-1) + A(3:2:Asize_x-1, 2*i-1))./4;
        p_new = ((-1)*(dati.c).*(dati.n).*(dati.dt)./(dati.length).*(Udelta) + p_old.*A_cap_old + A_cap_old-A_cap_new)./A_cap_new;
        p(1:Xsize-1,i) = p_new;
    %     p(Xsize, i) = p(Xsize-1, i);
        A_cap_old = A_cap_new;
        p_old = p_new;
        Pdelta = p(2:Xsize,i) - p(1:Xsize-1,i);
    
        
        figure(3);
        plot(x_dom,full(u(:,i)));
        title(['Volume velocity plot over time: ',num2str(t_dom(i)), ' s']);
        xlabel('Distance from the glottis [m]');
        ylabel('u [m^2]');
        ylim([-0.1,0.1]);
        hold on;
        plot(x_dom,sqrt(A(1:2:Asize_x, i*2-1)./pi),'r','Linewidth',2);
        plot(x_dom,-sqrt(A(1:2:Asize_x, i*2-1)./pi),'r','Linewidth',2);
        hold off
        pause(0.02);
    end
    
    Vout(vow,:) = u(Xsize, :);
    Vin= u(1,:);
    

    Y = fft(Vout(vow,:));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    Y = fft(Vin);
    P2 = abs(Y/L);
    P1_in = P2(1:L/2+1);
    P1_in(2:end-1) = 2*P1_in(2:end-1);

    H(vow, :)=P1./P1_in;
    H_log(vow, :)=20*log10(H(vow, :));    
    
    [~, locs] = findpeaks(H_log(vow, :));
    formants(vow,:) = f(locs(1:3));
    
end
    

figure('Name','Frequency Responce of vocal tract for differt vowels');
for vow=1:length(vows)
    plot(f,H_log(vow, :));
    xlim([0, 5e3]);
    hold on;
end
legend (['vowel ', vows(1)],['vowel ', vows(2)],['vowel ', vows(3)],['vowel ', vows(4)],['vowel ', vows(5)]);
ylabel('$20log_{10}(\left | H \right |)$', 'Interpreter','latex', 'FontSize',13);
xlabel('Frequency [Hz]');
hold off;

%%%%%%%%%%%%% formant deformation analysis %%%%%%%%%%%%%

% couples_formants = minDispVowelsPattern(formants);
% patternCost_formants = zeros(size(pattern,1),1);
% 
% for i = 1:size(pattern,1)
%     patternCost_formants(i,1) = couples_formants(pattern(i,1), pattern(i,2)) + couples_formants(pattern(i,2), pattern(i,3)) + couples_formants(pattern(i,3), pattern(i,4)) + couples_formants(pattern(i,4), pattern(i,5));
% end
% 
% minCost_formants = min(patternCost_formants(patternCost_formants>0), [], "all");
% minPattern_formants = find(patternCost_formants == minCost_formants);
% for i = 1:size(minPattern_formants,1)
%     i1 = pattern(minPattern_formants(i,1),1);
%     i2 = pattern(minPattern_formants(i,1),2);
%     i3 = pattern(minPattern_formants(i,1),3);
%     i4 = pattern(minPattern_formants(i,1),4);
%     i5 = pattern(minPattern_formants(i,1),5);
%     fprintf(['\n','min pattern for formants ', i, '= ',vows(i1), vows(i2), vows(i3), vows(i4), vows(i5) , '\n']);
% end

%%%%%%%%%%%%%%% Synthesis Sound %%%%%%%%%%%%%%%%%%%%%%%

for vow=1:5
    audio = Vout(vow,:);
    soundsc(audio, dati.fs);
end

end
