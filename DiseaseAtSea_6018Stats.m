load('2017001_6018_RRHmod_TV3.mat');

std(DiveStat.PDI/60,'omitnan')
std(DiveStat.PDI(DiveStat.SolarEl>=0)/60,'omitnan')
std(DiveStat.PDI(DiveStat.SolarEl<0)/60,'omitnan')

std(DiveStat.Dduration/60,'omitnan')
std(DiveStat.Dduration(DiveStat.SolarEl>=0)/60,'omitnan')
std(DiveStat.Dduration(DiveStat.SolarEl<0)/60,'omitnan')

std(DiveStat.Botttime/60,'omitnan')
std(DiveStat.Botttime(DiveStat.SolarEl>=0)/60,'omitnan')
std(DiveStat.Botttime(DiveStat.SolarEl<0)/60,'omitnan')

std(DiveStat.AscRate,'omitnan')
std(DiveStat.AscRate(DiveStat.SolarEl>=0),'omitnan')
std(DiveStat.AscRate(DiveStat.SolarEl<0),'omitnan')

std(DiveStat.DescRate,'omitnan')
std(DiveStat.DescRate(DiveStat.SolarEl>=0),'omitnan')
std(DiveStat.DescRate(DiveStat.SolarEl<0),'omitnan')

std(DiveStat.Maxdepth,'omitnan')
std(DiveStat.Maxdepth(DiveStat.SolarEl>=0),'omitnan')
std(DiveStat.Maxdepth(DiveStat.SolarEl<0),'omitnan')