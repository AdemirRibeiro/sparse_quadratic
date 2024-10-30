ηs = round.(output[:,5], sigdigits = 2);
ps = round.(output[:,6], sigdigits = 4);
ds = round.(output[:,7], sigdigits = 4);
ρs = round.(output[:,9], sigdigits = 2);
output[:,5] = ηs;
output[:,6] = ps;
output[:,7] = ds;
output[:,9] = ρs;

output_file = open("output_file.jl","w");
show(output_file, output);
close(output_file)
