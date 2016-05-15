i=1;
for SNR = 8:0.5:14
    M1 = 2;
    M2 = 2;
    test_dfe;
    P_GMDFE(i) = Pbit;
    P_theo(i) = Pbit_upper_bound;
    i=i+1;
end
i=1;
for SNR = 8:0.5:14
    M1 = 4;
    M2 = 0;
    test_dfe;
    P_GMLE(i) = Pbit;
    i=i+1;
end
i=1;
for SNR = 8:0.5:14
    M1 = 3;
    M2 = 4;
    test_aa_2;
    P_GMAA2(i) = Pbit;
    i=i+1;
end
i=1;
for SNR = 8:0.5:14
    M1 = 2;
    M2 = 4;
    test_aa;
    P_GMAA(i) = Pbit;
    i=i+1;
end

semilogy(linspace(8,14,13),P_GMAA2,'y');
hold on
semilogy(linspace(8,14,13),P_GMAA,'c');
hold on
semilogy(linspace(8,14,13),P_GMDFE,'r');
hold on
semilogy(linspace(8,14,13),P_GMLE,'g');
hold on
semilogy(linspace(8,14,13),P_theo,'b');

xlabel('SNR_(db)'); legend('AA + DFE @T/2','AA + DFE @T','MF + DFE','GM + LE','theoretical');
ylabel('Pbit') 