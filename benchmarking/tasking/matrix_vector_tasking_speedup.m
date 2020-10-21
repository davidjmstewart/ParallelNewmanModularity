clear all;
close all;

fig = figure;
matrix_sizes = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600];
T = [1: 16];

no_tasking_parallel_timings = [
0.000011490  0.000536640  0.001455860  0.000541430  0.001169260  0.002198460  0.004615750  0.006336860  0.001383060  0.000828970  0.001068510  0.001065800  0.001146610  0.001008810  0.001221530  0.001164140;
0.000044340  0.000469390  0.001301980  0.001791550  0.001108380  0.003220050  0.003312400  0.005804970  0.001473280  0.001026120  0.001014000  0.001279650  0.001135980  0.001244550  0.001225090  0.001272840;
0.000241030  0.000456040  0.001413750  0.000985280  0.001290850  0.003107240  0.006008550  0.006626480  0.001061600  0.000971220  0.000973880  0.001105430  0.001044730  0.001235750  0.001197340  0.001250690;
0.001343170  0.001046320  0.001358750  0.001210960  0.001287160  0.002185000  0.006143410  0.006370430  0.001992560  0.001577110  0.001413990  0.001372300  0.001422750  0.001542480  0.001520270  0.001731250;
0.002876740  0.002775130  0.002351010  0.002158230  0.003619600  0.003977430  0.006401420  0.006792260  0.002837820  0.002449940  0.003305850  0.002569540  0.002475820  0.002470120  0.002542930  0.002512490;
0.011150530  0.009065590  0.007647280  0.006933980  0.006950160  0.007191500  0.008981200  0.010741400  0.007284040  0.006940580  0.006568410  0.007037090  0.006970310  0.007393920  0.007072670  0.006830160;
0.051976710  0.032198570  0.026264610  0.024048170  0.024148720  0.024620590  0.027670040  0.027235360  0.024651460  0.024217940  0.023361080  0.024325570  0.025082320  0.024198240  0.024794290  0.025699790;
0.185167400  0.128539850  0.109714070  0.097075340  0.095108760  0.090413060  0.089130210  0.093387780  0.094776410  0.091498080  0.094078710  0.099926910  0.091120360  0.091702950  0.092504870  0.094084380;
0.763994250  0.491079500  0.433893020  0.411900190  0.367695040  0.355317330  0.351835830  0.350071120  0.372339510  0.356042010  0.349629690  0.354020290  0.349020070  0.353945040  0.362782520  0.356038130;]


no_tasking_sequential_averages = [
 0.000018840
 0.000093170
 0.000454820
 0.002528740
 0.006602870
 0.025474850
 0.114586420
 0.434647510
 1.610033300
]

with_tasking_parallel_timings = [
0.000014380  0.000162930  0.000009780  0.000675220  0.000029360  0.000688330  0.000007640  0.000661160  0.000006330  0.000912970  0.000011460  0.002407040  0.000387270  0.003460640  0.002547340  0.004142850;
0.000045090  0.000030480  0.000030180  0.000716550  0.000019010  0.000532360  0.000012320  0.001332570  0.000032550  0.001021710  0.000176160  0.003512170  0.000686210  0.005689980  0.003605980  0.005650630;
0.000231060  0.000185900  0.000156310  0.000540620  0.000107080  0.000529690  0.000098100  0.001374950  0.000458830  0.002779930  0.000727630  0.005880970  0.002995700  0.005878590  0.004440010  0.007259050;
0.001120990  0.001063320  0.001111330  0.000978890  0.001013600  0.001139000  0.000639740  0.001104000  0.000445140  0.002412420  0.000676640  0.003249440  0.000759410  0.005480490  0.003711430  0.006145680;
0.002484940  0.002673470  0.002406300  0.002454510  0.002098260  0.002053270  0.001828780  0.002107240  0.001425440  0.002613250  0.001426120  0.002827540  0.001370530  0.004109010  0.004316830  0.007278160;
0.011585630  0.011816540  0.011409240  0.008949490  0.008072220  0.006695170  0.006073270  0.006095180  0.005410730  0.006355100  0.005414350  0.006284400  0.005604930  0.009972210  0.009371450  0.010503680;
0.055721900  0.052722310  0.044655950  0.032721840  0.035657660  0.032073730  0.031196750  0.030379670  0.028444690  0.029183080  0.029420140  0.030393430  0.031384040  0.032905760  0.030570720  0.031072440;
0.174659400  0.196682930  0.168732360  0.131798400  0.142765350  0.113040120  0.104626790  0.085598530  0.094744220  0.095248460  0.097777480  0.094097820  0.096159910  0.092216520  0.085791930  0.088999590;
0.814605940  0.796531450  0.787852030  0.502393200  0.504733550  0.430740220  0.442288460  0.392531400  0.392045660  0.375313440  0.367180660  0.376853380  0.365792580  0.365175540  0.376519050  0.373828260;]


with_tasking_sequential_averages = [
0.000067140
0.000085690
0.000550550
0.001851920
0.006120390
0.025222180
0.127678340
0.441776610
1.665505360
]
for i = 1:size(no_tasking_parallel_timings, 1)
    subplot(5,2,i); 
    hold on;
    grid on;
%     xlim([1 16]);
%     ylim([0 6]);
    hold on;
    plot(T, no_tasking_parallel_timings(i,:).^-1 * no_tasking_sequential_averages(i), 'LineWidth', 2);
    plot(T, with_tasking_parallel_timings(i,:).^-1 * with_tasking_sequential_averages(i), 'LineWidth', 2);
    hold off
    title(sprintf('%d element Vector',matrix_sizes(i)));
end

handle=axes(fig,'visible','off'); 
handle.XLabel.Visible='on';
handle.YLabel.Visible='on';
handle.Title.Visible='on';
ylabel(handle,'Speedup');
xlabel(handle,'Number of threads');
title(handle,'Speedup of parallelising Matrix-Vector multiplication with SIMD support and tasks');
handle.Title.Position(2) = 1.05
grid on;

tic
VEC_SIZE = 500000000;
dot(rand(1,VEC_SIZE), rand(1,VEC_SIZE));
    toc