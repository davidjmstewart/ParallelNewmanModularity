clear all;
close all;

fig = figure;
matrix_sizes = [128, 256, 512, 1024, 2048, 4096];

parallel_timings_no_simd = [
    0.020185600  0.011102210  0.010713550  0.009270720  0.012251430  0.026384570  0.023694910  0.422775530;
    0.129443240  0.058825340  0.055001560  0.048942400  0.048647270  0.040422010  0.085471470  0.624313280;
    0.710216070  0.319322140  0.243795280  0.209328220  0.196237870  0.243558710  0.347685960  1.118053990;
    2.918158080  1.661382340  1.176851330  0.879510620  1.105529040  1.062844050  1.423657930  3.394205790;
    22.739670200  11.551353810  8.137235470  6.359376250  5.827440440  6.929228270  6.550827660  10.279084370;
    113.492606650  48.233910600  44.415683790  26.036232760  24.357439360  23.922651330  22.377171480  29.955789620  
];

parallel_timings_with_simd = [
    0.003967490  0.002452620  0.002334280  0.002835680  0.003167310  0.004779605  0.014810515  0.157470350;
    0.021284040  0.013979090  0.012096865  0.016438000  0.020345860  0.028354425  0.037038895  0.392328485;
    0.086523335  0.051200795  0.032284680  0.029365725  0.033995785  0.041194510  0.067170625  0.565975000;
    0.655887230  0.309871675  0.223601550  0.161743545  0.166930565  0.131436865  0.239247540  1.155852995;
    4.265759180  2.647362270  2.668531245  1.957684560  1.917737520  2.217248595  2.604632715  5.371673570;
    30.191197725  21.916348780  19.158561335  18.652140285  14.917684915  18.130009730  19.044290215  22.892084875 
];



T = [1: 8];

sequential_averages_simd = [
    0.007171645 
    0.060294735 
    0.333857580 
    1.616775470 
    9.286128890 
    72.838022155
];

sequential_averages_no_simd = [
    0.021040840 
    0.102657780 
    0.748604620 
    3.057870260 
    22.976747330 
    101.682787840 
];

T = [1: 8];


python_sequential_averages = [
0.027342939999999993, 0.12239112999999993, 0.49092758999999975, 1.6304126300000008, 6.983223840000005, 39.60936788999997];


for i = 1:size(parallel_timings_no_simd, 1)
    subplot(3,2,i); 
    hold on;
    grid on;
    grid minor;

    hold on;
    c_par = plot(T, parallel_timings_no_simd(i,:).^-1 * sequential_averages_no_simd(i), 'LineWidth', 2);
    p_par = plot(T, parallel_timings_no_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);
    c_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * sequential_averages_simd(i), 'LineWidth', 2);
    p_par_simd = plot(T, parallel_timings_with_simd(i,:).^-1 * python_sequential_averages(i), 'LineWidth', 2);

    hold off;

    title(sprintf('%d x %d matrix',matrix_sizes(i),matrix_sizes(i)),'interpreter','latex','fontsize',12);
    set(gca,'TickLabelInterpreter','latex')

end

hL = legend([c_par,p_par,c_par_simd,p_par_simd],{'Speedup (-O0 relative to sequential C)','Speedup (-O0 relative to Python)', 'Speedup with autovectorisation (-O3 relative to sequential C)', 'Speedup with autovectorisation (-O3 relative to Python)'},'interpreter','latex','fontsize',12);


handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.XLabel.Interpreter = 'latex';


handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup','interpreter','latex','fontsize',16);
xlabel(handle,'Number of threads','interpreter','latex','fontsize',16);
title(handle,'\textbf{Speedup of parallelising eigenpair calculations using autovectorisation \& OpenMP threads}','interpreter','latex','fontsize',12);
handle.Title.Position(2) = 1.05
