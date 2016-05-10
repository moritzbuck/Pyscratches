set terminal postscript color solid "Courier" 8
set output "bin18_reclust.ps"
set xtics rotate ( \
 "NODE_1_length_66333_cov_103.657_ID_1" 1, \
 "NODE_2_length_59894_cov_37.3764_ID_3" 66333, \
 "NODE_3_length_55612_cov_190.156_ID_5" 126226, \
 "NODE_4_length_53538_cov_33.1562_ID_7" 181837, \
 "NODE_5_length_49782_cov_29.9921_ID_9" 235374, \
 "NODE_6_length_48978_cov_42.2459_ID_11" 285155, \
 "NODE_7_length_45545_cov_78.8531_ID_13" 334132, \
 "NODE_8_length_41742_cov_30.7823_ID_15" 379676, \
 "NODE_9_length_40379_cov_89.0369_ID_17" 421417, \
 "NODE_45_length_17007_cov_103.352_ID_89" 461795, \
 "NODE_10_length_39219_cov_41.8163_ID_19" 478801, \
 "NODE_11_length_38499_cov_36.2448_ID_21" 518019, \
 "NODE_12_length_38432_cov_77.4041_ID_23" 556517, \
 "NODE_13_length_38421_cov_36.08_ID_25" 594948, \
 "NODE_14_length_36482_cov_56.1184_ID_27" 633368, \
 "NODE_15_length_36405_cov_56.4639_ID_29" 669849, \
 "NODE_128_length_6167_cov_158.599_ID_255" 706253, \
 "*NODE_198_length_3658_cov_123.974_ID_395" 712419, \
 "*NODE_173_length_4562_cov_140.176_ID_345" 716076, \
 "NODE_16_length_36302_cov_80.1348_ID_31" 720637, \
 "NODE_17_length_36124_cov_111.213_ID_33" 756938, \
 "NODE_18_length_35198_cov_38.7568_ID_35" 793061, \
 "NODE_97_length_7926_cov_32.4906_ID_193" 828258, \
 "NODE_175_length_4529_cov_115.071_ID_349" 836183, \
 "*NODE_180_length_4337_cov_63.4333_ID_359" 840711, \
 "NODE_49_length_15638_cov_42.5453_ID_97" 845047, \
 "*NODE_187_length_3997_cov_231.257_ID_373" 860684, \
 "*NODE_26_length_27074_cov_220.844_ID_51" 864680, \
 "NODE_226_length_2995_cov_45.4973_ID_451" 891754, \
 "NODE_19_length_34780_cov_20.6412_ID_37" 894748, \
 "NODE_20_length_34545_cov_28.4789_ID_39" 929527, \
 "NODE_167_length_4778_cov_20.8571_ID_333" 964071, \
 "NODE_21_length_34529_cov_54.5498_ID_41" 968847, \
 "NODE_185_length_4102_cov_121.119_ID_369" 1003375, \
 "NODE_22_length_33971_cov_69.7418_ID_43" 1007476, \
 "NODE_23_length_32531_cov_87.5592_ID_45" 1041446, \
 "NODE_76_length_9268_cov_64.4735_ID_151" 1073975, \
 "NODE_24_length_30529_cov_68.5953_ID_47" 1083242, \
 "*NODE_192_length_3840_cov_26.848_ID_383" 1113770, \
 "NODE_25_length_27226_cov_26.1728_ID_49" 1117609, \
 "NODE_82_length_8900_cov_30.4251_ID_163" 1144834, \
 "*NODE_135_length_5836_cov_29.6991_ID_269" 1153733, \
 "NODE_27_length_26814_cov_23.3162_ID_53" 1159568, \
 "*NODE_114_length_6770_cov_40.8034_ID_227" 1186381, \
 "NODE_28_length_26670_cov_37.228_ID_55" 1193150, \
 "*NODE_145_length_5509_cov_30.7406_ID_289" 1219819, \
 "NODE_112_length_6854_cov_21.9249_ID_223" 1225327, \
 "NODE_29_length_26205_cov_20.54_ID_57" 1232180, \
 "NODE_164_length_5005_cov_17.2045_ID_327" 1258384, \
 "NODE_30_length_25630_cov_36.9836_ID_59" 1263388, \
 "NODE_59_length_11521_cov_24.9368_ID_117" 1289017, \
 "NODE_190_length_3861_cov_57.3591_ID_379" 1300536, \
 "NODE_31_length_24909_cov_54.292_ID_61" 1304396, \
 "*NODE_38_length_18159_cov_62.2349_ID_75" 1329304, \
 "NODE_32_length_22506_cov_17.0397_ID_63" 1347462, \
 "*NODE_74_length_9420_cov_17.9695_ID_147" 1369967, \
 "NODE_33_length_21704_cov_13.6109_ID_65" 1379386, \
 "*NODE_58_length_11673_cov_264.416_ID_115" 1401089, \
 "NODE_34_length_21441_cov_263.165_ID_67" 1412763, \
 "NODE_35_length_20907_cov_95.3841_ID_69" 1434204, \
 "*NODE_174_length_4553_cov_106.475_ID_347" 1455112, \
 "NODE_155_length_5300_cov_113.365_ID_309" 1459664, \
 "NODE_36_length_18896_cov_20.4706_ID_71" 1464963, \
 "*NODE_131_length_6081_cov_40.9797_ID_261" 1483858, \
 "*NODE_218_length_3322_cov_35.9153_ID_435" 1489938, \
 "NODE_37_length_18462_cov_32.076_ID_73" 1493259, \
 "*NODE_182_length_4219_cov_35.4669_ID_363" 1511720, \
 "NODE_39_length_18002_cov_39.0522_ID_77" 1515937, \
 "NODE_40_length_17892_cov_33.7206_ID_79" 1533938, \
 "NODE_41_length_17407_cov_35.9377_ID_81" 1551829, \
 "NODE_42_length_17292_cov_54.152_ID_83" 1569232, \
 "*NODE_90_length_8491_cov_40.5853_ID_179" 1586523, \
 "NODE_43_length_17242_cov_80.4517_ID_85" 1595013, \
 "NODE_51_length_14630_cov_90.9911_ID_101" 1612255, \
 "NODE_91_length_8358_cov_132.522_ID_181" 1626884, \
 "NODE_44_length_17176_cov_120.552_ID_87" 1635241, \
 "NODE_46_length_16646_cov_15.3974_ID_91" 1652416, \
 "NODE_47_length_16074_cov_30.6975_ID_93" 1669061, \
 "NODE_80_length_9004_cov_29.3113_ID_159" 1685134, \
 "NODE_48_length_15697_cov_97.1941_ID_95" 1694137, \
 "*NODE_107_length_7248_cov_81.5595_ID_213" 1709833, \
 "*NODE_201_length_3587_cov_104.66_ID_401" 1717080, \
 "NODE_50_length_15170_cov_15.0753_ID_99" 1720666, \
 "*NODE_225_length_3012_cov_59.5302_ID_449" 1735835, \
 "NODE_84_length_8820_cov_33.3047_ID_167" 1738846, \
 "*NODE_195_length_3780_cov_39.9546_ID_389" 1747665, \
 "NODE_52_length_14250_cov_43.8339_ID_103" 1751444, \
 "NODE_53_length_13715_cov_30.9911_ID_105" 1765693, \
 "NODE_54_length_12778_cov_97.0846_ID_107" 1779408, \
 "NODE_55_length_12462_cov_17.7092_ID_109" 1792185, \
 "NODE_56_length_12066_cov_16.4283_ID_111" 1804646, \
 "NODE_57_length_11916_cov_12.1166_ID_113" 1816711, \
 "NODE_60_length_11283_cov_12.1063_ID_119" 1828626, \
 "NODE_61_length_10531_cov_35.4906_ID_121" 1839908, \
 "NODE_62_length_10530_cov_52.4296_ID_123" 1850438, \
 "NODE_63_length_10400_cov_35.6539_ID_125" 1860967, \
 "NODE_64_length_10321_cov_18.2854_ID_127" 1871366, \
 "NODE_65_length_10268_cov_84.4046_ID_129" 1881686, \
 "*NODE_143_length_5537_cov_77.7416_ID_285" 1891953, \
 "NODE_66_length_10192_cov_49.2487_ID_131" 1897489, \
 "NODE_132_length_5947_cov_86.5579_ID_263" 1907680, \
 "*NODE_188_length_3980_cov_69.6733_ID_375" 1913626, \
 "*NODE_236_length_2609_cov_128.324_ID_471" 1917605, \
 "NODE_67_length_10153_cov_19.0495_ID_133" 1920213, \
 "NODE_141_length_5652_cov_19.7681_ID_281" 1930365, \
 "NODE_68_length_9961_cov_22.4176_ID_135" 1936016, \
 "NODE_69_length_9829_cov_16.0126_ID_137" 1945976, \
 "NODE_70_length_9526_cov_37.8589_ID_139" 1955804, \
 "NODE_72_length_9514_cov_25.6831_ID_143" 1965329, \
 "NODE_71_length_9514_cov_13.2712_ID_141" 1974842, \
 "NODE_73_length_9442_cov_15.2898_ID_145" 1984355, \
 "NODE_75_length_9283_cov_110.661_ID_149" 1993796, \
 "NODE_77_length_9261_cov_15.797_ID_153" 2003078, \
 "NODE_78_length_9194_cov_26.9107_ID_155" 2012338, \
 "NODE_79_length_9164_cov_24.5417_ID_157" 2021531, \
 "*NODE_139_length_5669_cov_155.768_ID_277" 2030694, \
 "NODE_81_length_8911_cov_116.493_ID_161" 2036362, \
 "*NODE_207_length_3465_cov_117.915_ID_413" 2045272, \
 "NODE_83_length_8871_cov_14.3761_ID_165" 2048736, \
 "NODE_231_length_2716_cov_127.983_ID_461" 2057606, \
 "NODE_176_length_4525_cov_124.393_ID_351" 2060321, \
 "NODE_85_length_8767_cov_106.144_ID_169" 2064845, \
 "NODE_86_length_8733_cov_461.575_ID_171" 2073611, \
 "NODE_87_length_8733_cov_37.187_ID_173" 2082343, \
 "NODE_88_length_8724_cov_320.257_ID_175" 2091075, \
 "NODE_89_length_8570_cov_79.0291_ID_177" 2099798, \
 "NODE_92_length_8293_cov_28.4697_ID_183" 2108367, \
 "NODE_93_length_8178_cov_19.5497_ID_185" 2116659, \
 "NODE_94_length_8043_cov_108.736_ID_187" 2124835, \
 "NODE_95_length_8035_cov_22.2577_ID_189" 2132877, \
 "NODE_96_length_8024_cov_32.8105_ID_191" 2140911, \
 "NODE_203_length_3535_cov_160.314_ID_405" 2148934, \
 "NODE_98_length_7689_cov_180.765_ID_195" 2152468, \
 "*NODE_194_length_3813_cov_87.618_ID_387" 2160156, \
 "*NODE_212_length_3408_cov_92.6337_ID_423" 2163968, \
 "NODE_99_length_7631_cov_80.157_ID_197" 2167375, \
 "NODE_100_length_7582_cov_34.8692_ID_199" 2175005, \
 "NODE_101_length_7538_cov_22.1762_ID_201" 2182586, \
 "NODE_232_length_2712_cov_139.949_ID_463" 2190123, \
 "NODE_102_length_7536_cov_127.188_ID_203" 2192834, \
 "*NODE_178_length_4357_cov_183.488_ID_355" 2200369, \
 "NODE_103_length_7528_cov_48.9836_ID_205" 2204725, \
 "*NODE_117_length_6637_cov_57.8925_ID_233" 2212252, \
 "*NODE_230_length_2778_cov_82.8449_ID_459" 2218888, \
 "NODE_104_length_7478_cov_16.4625_ID_207" 2221666, \
 "NODE_105_length_7471_cov_31.1313_ID_209" 2229143, \
 "NODE_106_length_7449_cov_37.7897_ID_211" 2236613, \
 "NODE_205_length_3473_cov_45.4031_ID_409" 2244061, \
 "NODE_108_length_7191_cov_16.3642_ID_215" 2247533, \
 "NODE_109_length_7041_cov_66.6489_ID_217" 2254723, \
 "NODE_110_length_7016_cov_58.9265_ID_219" 2261763, \
 "*NODE_177_length_4494_cov_88.8454_ID_353" 2268778, \
 "NODE_111_length_6915_cov_20.6226_ID_221" 2273271, \
 "NODE_113_length_6777_cov_10.6255_ID_225" 2280185, \
 "NODE_115_length_6757_cov_108.899_ID_229" 2286961, \
 "*NODE_146_length_5503_cov_190.177_ID_291" 2293717, \
 "NODE_116_length_6754_cov_12.8173_ID_231" 2299219, \
 "NODE_118_length_6603_cov_34.0739_ID_235" 2305972, \
 "NODE_119_length_6599_cov_16.9494_ID_237" 2312574, \
 "NODE_120_length_6492_cov_18.0368_ID_239" 2319172, \
 "NODE_121_length_6465_cov_114.576_ID_241" 2325663, \
 "NODE_122_length_6399_cov_473.747_ID_243" 2332127, \
 "NODE_123_length_6394_cov_69.6687_ID_245" 2338525, \
 "NODE_124_length_6388_cov_22.9827_ID_247" 2344918, \
 "NODE_125_length_6375_cov_54.7704_ID_249" 2351305, \
 "NODE_126_length_6292_cov_42.0719_ID_251" 2357679, \
 "NODE_127_length_6248_cov_22.1731_ID_253" 2363970, \
 "NODE_129_length_6165_cov_29.488_ID_257" 2370217, \
 "NODE_140_length_5655_cov_144.58_ID_279" 2376381, \
 "NODE_130_length_6116_cov_139.402_ID_259" 2382035, \
 "NODE_133_length_5884_cov_12.4226_ID_265" 2388150, \
 "*NODE_161_length_5170_cov_101.083_ID_321" 2394032, \
 "NODE_134_length_5862_cov_322.135_ID_267" 2399201, \
 "NODE_137_length_5802_cov_21.2505_ID_273" 2405062, \
 "NODE_204_length_3476_cov_30.9217_ID_407" 2410863, \
 "*NODE_202_length_3540_cov_14.4811_ID_403" 2414338, \
 "NODE_136_length_5802_cov_11.8449_ID_271" 2417877, \
 "NODE_138_length_5704_cov_75.3494_ID_275" 2423678, \
 "NODE_142_length_5647_cov_90.1104_ID_283" 2429381, \
 "NODE_144_length_5516_cov_20.4041_ID_287" 2435027, \
 "NODE_147_length_5455_cov_17.6359_ID_293" 2440542, \
 "NODE_148_length_5441_cov_115.461_ID_295" 2445996, \
 "NODE_149_length_5430_cov_69.2046_ID_297" 2451436, \
 "NODE_150_length_5430_cov_19.2146_ID_299" 2456865, \
 "NODE_151_length_5365_cov_42.6785_ID_301" 2462294, \
 "NODE_152_length_5352_cov_352.825_ID_303" 2467658, \
 "NODE_153_length_5312_cov_64.51_ID_305" 2473009, \
 "NODE_154_length_5305_cov_12.5197_ID_307" 2478320, \
 "NODE_156_length_5280_cov_13.3982_ID_311" 2483624, \
 "NODE_157_length_5261_cov_18.2591_ID_313" 2488903, \
 "NODE_158_length_5246_cov_24.509_ID_315" 2494163, \
 "NODE_159_length_5236_cov_61.8228_ID_317" 2499408, \
 "*NODE_211_length_3413_cov_54.6628_ID_421" 2504643, \
 "NODE_193_length_3815_cov_81.9235_ID_385" 2508055, \
 "NODE_160_length_5221_cov_329.894_ID_319" 2511869, \
 "NODE_162_length_5068_cov_41.449_ID_323" 2517089, \
 "NODE_163_length_5024_cov_14.6517_ID_325" 2522156, \
 "NODE_165_length_4927_cov_23.8414_ID_329" 2527179, \
 "NODE_166_length_4830_cov_56.0903_ID_331" 2532105, \
 "NODE_168_length_4718_cov_167.852_ID_335" 2536934, \
 "*NODE_196_length_3775_cov_141.002_ID_391" 2541651, \
 "*NODE_222_length_3183_cov_65.735_ID_443" 2545425, \
 "NODE_169_length_4707_cov_41.6739_ID_337" 2548607, \
 "NODE_170_length_4669_cov_86.404_ID_339" 2553313, \
 "NODE_171_length_4590_cov_51.9597_ID_341" 2557981, \
 "NODE_172_length_4585_cov_194.417_ID_343" 2562570, \
 "NODE_179_length_4352_cov_85.5247_ID_357" 2567154, \
 "NODE_181_length_4259_cov_143.278_ID_361" 2571505, \
 "NODE_191_length_3846_cov_156.191_ID_381" 2575763, \
 "NODE_183_length_4187_cov_170.924_ID_365" 2579608, \
 "NODE_184_length_4106_cov_30.0655_ID_367" 2583794, \
 "NODE_186_length_4030_cov_799.18_ID_371" 2587899, \
 "NODE_215_length_3357_cov_439.734_ID_429" 2591928, \
 "NODE_189_length_3895_cov_31.5236_ID_377" 2595284, \
 "*NODE_214_length_3366_cov_146.379_ID_427" 2599178, \
 "NODE_197_length_3749_cov_188.86_ID_393" 2602543, \
 "NODE_199_length_3653_cov_166.87_ID_397" 2606291, \
 "NODE_200_length_3629_cov_69.6357_ID_399" 2609943, \
 "NODE_206_length_3470_cov_135.036_ID_411" 2613571, \
 "NODE_208_length_3455_cov_730.621_ID_415" 2617040, \
 "NODE_223_length_3055_cov_194.051_ID_445" 2620494, \
 "NODE_209_length_3431_cov_139.031_ID_417" 2623548, \
 "NODE_210_length_3424_cov_50.2931_ID_419" 2626978, \
 "NODE_213_length_3388_cov_547.21_ID_425" 2630401, \
 "NODE_216_length_3347_cov_231.909_ID_431" 2633788, \
 "NODE_219_length_3305_cov_413.311_ID_437" 2637134, \
 "NODE_220_length_3232_cov_18.0482_ID_439" 2640438, \
 "NODE_221_length_3221_cov_23.6027_ID_441" 2643666, \
 "*NODE_235_length_2652_cov_77.1973_ID_469" 2646886, \
 "NODE_224_length_3018_cov_42.0296_ID_447" 2649537, \
 "NODE_227_length_2971_cov_37.4229_ID_453" 2652554, \
 "NODE_228_length_2936_cov_323.836_ID_455" 2655524, \
 "NODE_229_length_2883_cov_47.4629_ID_457" 2658459, \
 "NODE_233_length_2683_cov_21.5944_ID_465" 2661341, \
 "NODE_234_length_2679_cov_135.603_ID_467" 2664023, \
 "NODE_217_length_3331_cov_162.805_ID_433" 2666701, \
 "" 2670266 \
)
set ytics ( \
 "contig-319" 1, \
 "contig-2000391" 56154, \
 "*contig-254" 66298, \
 "*contig-596" 126148, \
 "contig-14000434" 185501, \
 "*contig-1000092" 238724, \
 "contig-1000034" 259213, \
 "*contig-3000727" 277639, \
 "*contig-77" 291355, \
 "contig-109" 341054, \
 "*contig-2000582" 388302, \
 "contig-1000331" 393862, \
 "*contig-16000658" 411951, \
 "contig-111" 431056, \
 "contig-18000132" 439763, \
 "contig-4000064" 488927, \
 "*contig-2000640" 528107, \
 "*contig-3000813" 566561, \
 "contig-174" 604944, \
 "*contig-13000746" 643347, \
 "*contig-1000012" 680403, \
 "contig-2000476" 714508, \
 "contig-235" 749623, \
 "*contig-662" 786698, \
 "contig-989" 824057, \
 "contig-5000090" 859646, \
 "*contig-175" 891099, \
 "*contig-295" 942054, \
 "contig-742" 961894, \
 "*contig-12000232" 977885, \
 "contig-794" 985004, \
 "contig-15000366" 1018779, \
 "*contig-1000048" 1042294, \
 "contig-916" 1054625, \
 "contig-584" 1067005, \
 "contig-470" 1101351, \
 "*contig-12000833" 1133598, \
 "*contig-5000265" 1141964, \
 "*contig-504" 1149942, \
 "*contig-9000901" 1186396, \
 "contig-21000286" 1225679, \
 "*contig-1000706" 1246421, \
 "contig-15000317" 1254107, \
 "contig-1005" 1269338, \
 "*contig-15000728" 1282775, \
 "*contig-726" 1291833, \
 "contig-507" 1297364, \
 "*contig-4000848" 1330563, \
 "*contig-903" 1338497, \
 "*contig-538" 1376964, \
 "*contig-13000227" 1413960, \
 "contig-17000594" 1462889, \
 "*contig-493" 1483564, \
 "*contig-763" 1495480, \
 "*contig-446" 1517148, \
 "contig-16000085" 1552219, \
 "*contig-1000651" 1566253, \
 "contig-483" 1579547, \
 "*contig-16000062" 1604794, \
 "*contig-17000898" 1625837, \
 "*contig-766" 1651798, \
 "contig-511" 1658384, \
 "contig-1000925" 1676368, \
 "contig-2000780" 1683341, \
 "*contig-249" 1689133, \
 "contig-404" 1694359, \
 "*contig-4000999" 1711741, \
 "contig-12000394" 1736328, \
 "*contig-14000732" 1747317, \
 "*contig-11000312" 1767049, \
 "contig-9000872" 1778677, \
 "contig-12000096" 1784847, \
 "contig-163" 1791039, \
 "*contig-896" 1810723, \
 "*contig-17000086" 1845945, \
 "*contig-13000089" 1862230, \
 "*contig-14000790" 1889883, \
 "*contig-6000443" 1913478, \
 "*contig-2000651" 1935065, \
 "contig-14000406" 1956997, \
 "contig-11000465" 1963444, \
 "contig-1000845" 1972706, \
 "contig-8000580" 1978642, \
 "*contig-3000161" 2009569, \
 "*contig-14000599" 2022322, \
 "contig-12000014" 2029955, \
 "contig-8000274" 2057699, \
 "*contig-721" 2070185, \
 "*contig-2000357" 2075386, \
 "*contig-842" 2082270, \
 "contig-0" 2088922, \
 "contig-14000127" 2094159, \
 "contig-3000673" 2103954, \
 "*contig-14000546" 2129436, \
 "*contig-13000761" 2140177, \
 "contig-545" 2151143, \
 "*contig-5000598" 2161425, \
 "*contig-2000705" 2168471, \
 "*contig-18000675" 2189846, \
 "*contig-10000471" 2206327, \
 "*contig-17000958" 2214272, \
 "*contig-783" 2219892, \
 "contig-12000503" 2230551, \
 "*contig-113" 2245975, \
 "contig-3000147" 2257009, \
 "*contig-12000986" 2266494, \
 "contig-15000177" 2275964, \
 "*contig-8000577" 2285436, \
 "*contig-2000756" 2294854, \
 "*contig-1000410" 2309021, \
 "*contig-951" 2318248, \
 "contig-6000727" 2327614, \
 "contig-6000842" 2336755, \
 "contig-527" 2382901, \
 "*contig-1000650" 2391575, \
 "contig-1000883" 2410933, \
 "contig-338" 2419908, \
 "*contig-11000533" 2428601, \
 "contig-830" 2437277, \
 "contig-7000488" 2445771, \
 "contig-1000799" 2453926, \
 "*contig-11000098" 2462067, \
 "contig-603" 2474432, \
 "*contig-14000823" 2482487, \
 "*contig-495" 2490478, \
 "*contig-1000146" 2523543, \
 "contig-646" 2545121, \
 "*contig-458" 2552679, \
 "*contig-92" 2560759, \
 "contig-302" 2581583, \
 "*contig-1000697" 2588010, \
 "*contig-2000270" 2604711, \
 "*contig-5000904" 2612179, \
 "contig-658" 2619605, \
 "*contig-1000205" 2634394, \
 "contig-2000247" 2641544, \
 "*contig-10000767" 2649201, \
 "*contig-1000305" 2661957, \
 "contig-436" 2668828, \
 "contig-1000287" 2675614, \
 "*contig-2000209" 2681063, \
 "contig-2000102" 2710371, \
 "*contig-8000613" 2716966, \
 "*contig-5000972" 2723694, \
 "*contig-6000343" 2730300, \
 "contig-11000589" 2736630, \
 "*contig-1000454" 2743948, \
 "contig-2000445" 2750389, \
 "*contig-7000756" 2756758, \
 "*contig-1000211" 2763396, \
 "contig-3000454" 2769803, \
 "contig-1000328" 2776140, \
 "contig-1000818" 2782296, \
 "contig-1000367" 2788541, \
 "contig-6000473" 2794681, \
 "contig-8000008" 2822599, \
 "contig-274" 2850523, \
 "*contig-12000241" 2856568, \
 "*contig-5001008" 2885930, \
 "contig-1000570" 2892053, \
 "*contig-15000630" 2897672, \
 "*contig-4000408" 2903055, \
 "*contig-124" 2911429, \
 "contig-1013" 2917237, \
 "*contig-12000928" 2922893, \
 "*contig-149" 2930753, \
 "contig-7000074" 2936480, \
 "*contig-667" 2941912, \
 "contig-8000727" 2947630, \
 "*contig-12000998" 2953011, \
 "contig-9000350" 2963274, \
 "*contig-628" 2968689, \
 "*contig-12000975" 2974895, \
 "*contig-146" 2980198, \
 "*contig-1000411" 2985274, \
 "*contig-13000492" 2990577, \
 "*contig-1000500" 2995823, \
 "contig-8000897" 3001085, \
 "*contig-8001015" 3008254, \
 "contig-11000317" 3020908, \
 "contig-2000891" 3029427, \
 "*contig-5000315" 3034452, \
 "*contig-981" 3039453, \
 "*contig-11000516" 3044943, \
 "*contig-710" 3055925, \
 "contig-3000791" 3067845, \
 "*contig-4000363" 3078793, \
 "*contig-3000485" 3102114, \
 "contig-6000273" 3108488, \
 "contig-1000824" 3129504, \
 "contig-2000153" 3142696, \
 "contig-894" 3151587, \
 "*contig-1000182" 3157155, \
 "contig-2000288" 3190669, \
 "contig-58" 3228811, \
 "contig-13000046" 3234472, \
 "contig-2000519" 3249101, \
 "*contig-12000949" 3257913, \
 "*contig-2000097" 3269938, \
 "*contig-12000954" 3290296, \
 "*contig-10000264" 3302869, \
 "contig-10000328" 3329465, \
 "contig-6000675" 3346501, \
 "*contig-3000923" 3356787, \
 "contig-1000931" 3367176, \
 "contig-9000547" 3373328, \
 "contig-1000638" 3407351, \
 "*contig-261" 3413836, \
 "contig-1000806" 3420198, \
 "contig-16000594" 3426046, \
 "contig-6000760" 3433138, \
 "*contig-865" 3441212, \
 "*contig-15000260" 3446340, \
 "contig-5000396" 3460745, \
 "*contig-8000535" 3465797, \
 "contig-2000034" 3474114, \
 "*contig-5000803" 3479260, \
 "contig-1000151" 3485341, \
 "contig-1000874" 3491314, \
 "contig-11000455" 3500547, \
 "contig-10000556" 3508312, \
 "contig-2000807" 3515828, \
 "contig-3000492" 3521883, \
 "contig-2000906" 3529158, \
 "contig-11000895" 3536375, \
 "contig-2000608" 3541918, \
 "contig-7000040" 3551678, \
 "contig-2000864" 3558347, \
 "contig-2000833" 3563939, \
 "contig-1000856" 3575874, \
 "contig-8000525" 3581632, \
 "contig-3000178" 3587376, \
 "contig-7000650" 3601878, \
 "contig-21000556" 3617221, \
 "contig-6000735" 3622380, \
 "" 3629159 \
)
set size 3,3
set grid
unset key
set border 0
set tics scale 0
set xlabel "REF"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:2670266]
set yrange [1:3629159]
set zrange [0:100]
set colorbox default
set cblabel "%similarity"
set cbrange [0:100]
set cbtics 20
set pm3d map
set palette model RGB defined ( \
  0 "#000000", \
  4 "#DD00DD", \
  6 "#0000DD", \
  7 "#00DDDD", \
  8 "#00DD00", \
  9 "#DDDD00", \
 10 "#DD0000"  \
)
set style line 1  palette lw 4 pt 6 ps 0.5
set style line 2  palette lw 4 pt 6 ps 0.5
set style line 3  palette lw 4 pt 6 ps 0.5
splot \
 "bin18_reclust.fplot" title "FWD" w l ls 1, \
 "bin18_reclust.rplot" title "REV" w l ls 2
