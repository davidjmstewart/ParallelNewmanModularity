clear all;
close all;

fig = figure;
matrix_sizes = [128, 256, 512, 1024, 2048, 4096, 8192, 16384];
% matrix_sizes = [1000000, 2000000, 3000000, 4000000];

parallel_timings_no_simd = [0.000049360  0.000713650  0.000631820  0.000916120  0.000890430  0.001971290  0.003043060  0.007647590 
0.000128760  0.001100380  0.001231370  0.002039780  0.002096400  0.004303130  0.007214820  0.010544880 
0.000718460  0.000687270  0.000909790  0.001142110  0.001436060  0.003223050  0.005432640  0.009145190 
0.002617670  0.001887430  0.001778370  0.001671050  0.002790780  0.003247430  0.003616580  0.007271120 
0.010399030  0.007451670  0.007066100  0.005850140  0.005841580  0.006304490  0.007835250  0.010624150 
0.044106600  0.029603090  0.023806290  0.020113950  0.019077340  0.019816350  0.019446320  0.022429880 
0.168941910  0.105218830  0.091166860  0.084881690  0.079223090  0.077749560  0.081812080  0.079671270 
0.663622500  0.412277470  0.352370110  0.326236490  0.316735790  0.311098190  0.302963940  0.310808110  ];

parallel_timings_with_simd = [ 0.000024310  0.000874800  0.000898490  0.001174240  0.001565370  0.002542570  0.004728760  0.006606670
 0.000123100  0.001780800  0.001331150  0.002258110  0.003416100  0.005117890  0.008763020  0.011863240
 0.000716830  0.000667460  0.001109050  0.001187780  0.002177000  0.003709000  0.004419020  0.007594160
 0.002551040  0.002074410  0.002464100  0.003047960  0.002268910  0.004517550  0.005842690  0.007953620
 0.008516630  0.007435820  0.007333710  0.006626040  0.006620370  0.006730940  0.008673720  0.010526490
 0.031229470  0.019279220  0.016812180  0.016021710  0.018891560  0.020317710  0.019497570  0.021395680
 0.121268680  0.090828280  0.076048850  0.077145530  0.071672710  0.075138120  0.072249590  0.077476640
 0.527066970  0.352656350  0.325065800  0.322691190  0.313213060  0.321134100  0.305272380  0.312454850 ];



T = [1: 8];

sequential_averages_no_simd = [
0.000033590 
0.000374240 
0.000984970 
0.003217040 
0.011640270 
0.042096940 
0.172534180 
0.680300710 
];


T = [1: 8];


for i = 1:size(parallel_timings_no_simd, 1)
    subplot(4,2,i); 
    hold on;
    grid on;
grid minor;
%     xlim([1 16]);
%     ylim([0 6]);
    hold on;
    c_par = plot(T, parallel_timings_with_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
    p_par = plot(T, parallel_timings_no_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
%     c_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
%     p_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);

    hold off;
%     set(gca, 'YScale', 'log')

    title(sprintf('%d element Vector',matrix_sizes(i)),'interpreter','latex','fontsize',12);
    set(gca,'TickLabelInterpreter','latex')

end

% hL = legend([c_par,p_par,c_par_simd,p_par_simd],{'Speedup (relative to sequential C)','Speedup (relative to Python)', 'Speedup with autovectorisation (relative to sequential C)', 'Speedup with autovectorisation (relative to Python)'},'interpreter','latex','fontsize',12);
handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.XLabel.Interpreter = 'latex';


handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup','interpreter','latex','fontsize',16);
xlabel(handle,'Number of threads','interpreter','latex','fontsize',16);
title(handle,'\textbf{Speedup of parallelising eigenpair calculations using autovectorisation \& OpenMP threads}','interpreter','latex','fontsize',12);
handle.Title.Position(2) = 1.05
