clear all;
close all;

fig = figure;
matrix_sizes = [128, 256, 512, 1024, 2048, 4096];
% matrix_sizes = [1000000, 2000000, 3000000, 4000000];

parallel_timings_no_simd = [0.002234070  0.000860010  0.001023785  0.001195810  0.001831060  0.003917260  0.010468015  0.060036250;
0.004865040  0.002016685  0.002342600  0.001558590  0.002185005  0.004756740  0.009758295  0.055143235;
0.013415185  0.006359940  0.004673390  0.003699355  0.003579310  0.003634385  0.004582165  0.051254820;
0.049079030  0.025823980  0.018595900  0.014118725  0.013079430  0.013276385  0.017795035  0.073375830;
0.174042270  0.091002175  0.064370630  0.050796330  0.044126725  0.038856430  0.051004800  0.097289785;
0.679244900  0.340198245  0.234881285  0.186297140  0.155703940  0.142103060  0.141232880  0.193348295; ];

parallel_timings_with_simd = [0.001363080  0.000843090  0.000565800  0.000756500  0.000991610  0.002179445  0.006567545  0.044328205;
0.001523360  0.000713600  0.001127765  0.001302130  0.001956140  0.002583240  0.007897505  0.060372835;
0.003151560  0.001233070  0.001407440  0.001348085  0.002184085  0.004303475  0.010724210  0.046791490;
0.010372015  0.006246125  0.005751550  0.004448820  0.003609430  0.003633125  0.006114680  0.045151040;
0.037017210  0.026520715  0.021984310  0.020463515  0.022294250  0.018845670  0.023805090  0.064307750;
0.129931610  0.089414740  0.076108410  0.073708230  0.069604460  0.068723160  0.070520780  0.118061550; ];



T = [1: 8];

sequential_averages_no_simd = [
0.002294190
0.004500290
0.013623485
0.049738380
0.174274960
0.673365850
];


T = [1: 8];


python_sequential_averages = [
0.027342939999999993, 0.12239112999999993, 0.49092758999999975, 1.6304126300000008, 6.983223840000005, 39.60936788999997];


for i = 1:size(parallel_timings_no_simd, 1)
    subplot(3,2,i); 
    hold on;
    grid on;
grid minor;
%     xlim([1 16]);
%     ylim([0 6]);
    hold on;
    c_par = plot(T, parallel_timings_no_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
    p_par = plot(T, parallel_timings_no_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);
    c_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
    p_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);

    hold off;
    set(gca, 'YScale', 'log')

    title(sprintf('%d element Vector',matrix_sizes(i)),'interpreter','latex','fontsize',12);
    set(gca,'TickLabelInterpreter','latex')

end

hL = legend([c_par,p_par,c_par_simd,p_par_simd],{'Speedup (relative to sequential C)','Speedup (relative to Python)', 'Speedup with autovectorisation (relative to sequential C)', 'Speedup with autovectorisation (relative to Python)'},'interpreter','latex','fontsize',12);
handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.XLabel.Interpreter = 'latex';


handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup','interpreter','latex','fontsize',16);
xlabel(handle,'Number of threads','interpreter','latex','fontsize',16);
title(handle,'\textbf{Speedup of parallelising eigenpair calculations using autovectorisation \& OpenMP threads}','interpreter','latex','fontsize',12);
handle.Title.Position(2) = 1.05
