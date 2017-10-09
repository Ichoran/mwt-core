# Copyright Nicholas A. Swierczek, Rex A. Kerr and HHMI Janelia, 2007-2015.
# Copyright Rex A. Kerr and Calico Life Sciences, 2015-2016.
# This file is covered by the LGPL 2.1 license.

SHELL=/bin/bash

FLAGS = -Wall -O2 -fno-strict-aliasing -ggdb3

ifeq ($(OS),Windows_NT)
    CC = mingw32-g++ -std=gnu++11
    OUTDIR = lib
    TGT = -DWINDOWS -DENABLE_SIMD
else
    CC = g++ -std=gnu++11
#   CC = clang++ -std=c++11 -DCLANG_WORKAROUND
    TGT = -DLINUX -DENABLE_SIMD
endif

UNIT = -DUNIT_TEST_OWNER
all: unit_geometry unit_lists unit_storage unit_image unit_align unit_blob unit_model unit_library

unit_geometry: makefile MWT_Geometry.h MWT_Geometry.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_geometry MWT_Geometry.cc

MWT_Geometry.o: makefile MWT_Geometry.h MWT_Geometry.cc
	$(CC) $(FLAGS) $(TGT) -o MWT_Geometry.o MWT_Geometry.cc
	
unit_lists: makefile MWT_Lists.h MWT_Lists.cc MWT_Geometry.h
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_lists MWT_Lists.cc

MWT_Lists.o: makefile MWT_Lists.h MWT_Lists.cc MWT_Geometry.h
	$(CC) $(FLAGS) $(TGT) -o MWT_Lists.o MWT_Lists.cc
	
unit_storage: makefile MWT_Storage.h MWT_Storage.cc MWT_Lists.h
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_storage MWT_Storage.cc

MWT_Image.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Image.o MWT_Image.cc

unit_image: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_image MWT_Image.cc

MWT_Align.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Align.h MWT_Align.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Align.o MWT_Align.cc

unit_align: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Align.h MWT_Align.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_align MWT_Align.cc MWT_Image.o

MWT_Blob.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Align.h MWT_Blob.h MWT_Blob.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Blob.o MWT_Blob.cc

unit_blob: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Align.o MWT_Blob.h MWT_Blob.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_blob MWT_Blob.cc MWT_Image.o
	
MWT_Model.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Model.h MWT_Model.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Model.o MWT_Model.cc
	
unit_model: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Model.h MWT_Model.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_model MWT_Model.cc MWT_Image.o

unit_library: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Align.h MWT_Align.o MWT_Blob.h MWT_Blob.o MWT_Model.h MWT_Model.o MWT_Library.h MWT_Library.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_library MWT_Library.cc MWT_Image.o MWT_Align.o MWT_Blob.o MWT_Model.o

MWT_Library.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Align.h MWT_Blob.h MWT_Model.h MWT_Library.h MWT_Library.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Library.o MWT_Library.cc

mwt_bench: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Align.h MWT_Blob.h MWT_Library.h MWT_Model.h MWT_Bench.cc MWT_Library.o MWT_Model.o MWT_Blob.o MWT_Image.o
	$(CC) $(FLAGS) $(TGT) -o mwt_bench MWT_Bench.cc MWT_Library.o MWT_Model.o MWT_Blob.o MWT_Align.o MWT_Image.o

ifeq ($(OS),Windows_NT)
    DLL: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Align.h MWT_Align.o MWT_Blob.h MWT_Blob.o MWT_Model.h MWT_Model.o MWT_Library.h MWT_Library.o MWT_DLL.h MWT_DLL.cc
		mkdir -p lib; $(CC) -shared -static -static-libstdc++ -static-libgcc $(FLAGS) $(TGT) -o $(OUTDIR)/MWT.dll MWT_DLL.cc MWT_Image.o MWT_Align.o MWT_Blob.o MWT_Model.o MWT_Library.o
endif

test: all
	./unit_geometry
	./unit_lists
	./unit_storage
	./unit_image
	./unit_align
	./unit_blob
	./unit_model
	./unit_library

bench: mwt_bench
	time ./mwt_bench fake-data/000.raw fake-data/001.raw fake-data/002.raw fake-data/003.raw fake-data/004.raw fake-data/005.raw fake-data/006.raw fake-data/007.raw fake-data/008.raw fake-data/009.raw fake-data/010.raw fake-data/011.raw fake-data/012.raw fake-data/013.raw fake-data/014.raw fake-data/015.raw fake-data/016.raw fake-data/017.raw fake-data/018.raw fake-data/019.raw fake-data/020.raw fake-data/021.raw fake-data/022.raw fake-data/023.raw fake-data/024.raw fake-data/025.raw fake-data/026.raw fake-data/027.raw fake-data/028.raw fake-data/029.raw fake-data/030.raw fake-data/031.raw fake-data/032.raw fake-data/033.raw fake-data/034.raw fake-data/035.raw fake-data/036.raw fake-data/037.raw fake-data/038.raw fake-data/039.raw fake-data/040.raw fake-data/041.raw fake-data/042.raw fake-data/043.raw fake-data/044.raw fake-data/045.raw fake-data/046.raw fake-data/047.raw fake-data/048.raw fake-data/049.raw fake-data/050.raw fake-data/051.raw fake-data/052.raw fake-data/053.raw fake-data/054.raw fake-data/055.raw fake-data/056.raw fake-data/057.raw fake-data/058.raw fake-data/059.raw fake-data/060.raw fake-data/061.raw fake-data/062.raw fake-data/063.raw fake-data/064.raw fake-data/065.raw fake-data/066.raw fake-data/067.raw fake-data/068.raw fake-data/069.raw fake-data/070.raw fake-data/071.raw fake-data/072.raw fake-data/073.raw fake-data/074.raw fake-data/075.raw fake-data/076.raw fake-data/077.raw fake-data/078.raw fake-data/079.raw fake-data/080.raw fake-data/081.raw fake-data/082.raw fake-data/083.raw fake-data/084.raw fake-data/085.raw fake-data/086.raw fake-data/087.raw fake-data/088.raw fake-data/089.raw fake-data/090.raw fake-data/091.raw fake-data/092.raw fake-data/093.raw fake-data/094.raw fake-data/095.raw fake-data/096.raw fake-data/097.raw fake-data/098.raw fake-data/099.raw fake-data/100.raw fake-data/101.raw fake-data/102.raw fake-data/103.raw fake-data/104.raw fake-data/105.raw fake-data/106.raw fake-data/107.raw fake-data/108.raw fake-data/109.raw fake-data/110.raw fake-data/111.raw fake-data/112.raw fake-data/113.raw fake-data/114.raw fake-data/115.raw fake-data/116.raw fake-data/117.raw fake-data/118.raw fake-data/119.raw fake-data/120.raw fake-data/121.raw fake-data/122.raw fake-data/123.raw fake-data/124.raw fake-data/125.raw fake-data/126.raw fake-data/127.raw fake-data/128.raw fake-data/129.raw fake-data/130.raw fake-data/131.raw fake-data/132.raw fake-data/133.raw fake-data/134.raw fake-data/135.raw fake-data/136.raw fake-data/137.raw fake-data/138.raw fake-data/139.raw fake-data/140.raw fake-data/141.raw fake-data/142.raw fake-data/143.raw fake-data/144.raw fake-data/145.raw fake-data/146.raw fake-data/147.raw fake-data/148.raw fake-data/149.raw fake-data/150.raw fake-data/151.raw fake-data/152.raw fake-data/153.raw fake-data/154.raw fake-data/155.raw fake-data/156.raw fake-data/157.raw fake-data/158.raw fake-data/159.raw fake-data/160.raw fake-data/161.raw fake-data/162.raw fake-data/163.raw fake-data/164.raw fake-data/165.raw fake-data/166.raw fake-data/167.raw fake-data/168.raw fake-data/169.raw fake-data/170.raw fake-data/171.raw fake-data/172.raw fake-data/173.raw fake-data/174.raw fake-data/175.raw fake-data/176.raw fake-data/177.raw fake-data/178.raw fake-data/179.raw fake-data/180.raw fake-data/181.raw fake-data/182.raw fake-data/183.raw fake-data/184.raw fake-data/185.raw fake-data/186.raw fake-data/187.raw fake-data/188.raw fake-data/189.raw fake-data/190.raw fake-data/191.raw fake-data/192.raw fake-data/193.raw fake-data/194.raw fake-data/195.raw fake-data/196.raw fake-data/197.raw fake-data/198.raw fake-data/199.raw fake-data/200.raw fake-data/201.raw fake-data/202.raw fake-data/203.raw fake-data/204.raw fake-data/205.raw fake-data/206.raw fake-data/207.raw fake-data/208.raw fake-data/209.raw fake-data/210.raw fake-data/211.raw fake-data/212.raw fake-data/213.raw fake-data/214.raw fake-data/215.raw fake-data/216.raw fake-data/217.raw fake-data/218.raw fake-data/219.raw fake-data/220.raw fake-data/221.raw fake-data/222.raw fake-data/223.raw fake-data/224.raw fake-data/225.raw fake-data/226.raw fake-data/227.raw fake-data/228.raw fake-data/229.raw fake-data/230.raw fake-data/231.raw fake-data/232.raw fake-data/233.raw fake-data/234.raw fake-data/235.raw fake-data/236.raw fake-data/237.raw fake-data/238.raw fake-data/239.raw fake-data/240.raw fake-data/241.raw fake-data/242.raw fake-data/243.raw fake-data/244.raw fake-data/245.raw fake-data/246.raw fake-data/247.raw fake-data/248.raw fake-data/249.raw fake-data/250.raw fake-data/251.raw fake-data/252.raw fake-data/253.raw fake-data/254.raw fake-data/255.raw fake-data/256.raw fake-data/257.raw fake-data/258.raw fake-data/259.raw fake-data/260.raw fake-data/261.raw fake-data/262.raw fake-data/263.raw fake-data/264.raw fake-data/265.raw fake-data/266.raw fake-data/267.raw fake-data/268.raw fake-data/269.raw fake-data/270.raw fake-data/271.raw fake-data/272.raw fake-data/273.raw fake-data/274.raw fake-data/275.raw fake-data/276.raw fake-data/277.raw fake-data/278.raw fake-data/279.raw fake-data/280.raw fake-data/281.raw fake-data/282.raw fake-data/283.raw fake-data/284.raw fake-data/285.raw fake-data/286.raw fake-data/287.raw fake-data/288.raw fake-data/289.raw fake-data/290.raw fake-data/291.raw fake-data/292.raw fake-data/293.raw fake-data/294.raw fake-data/295.raw fake-data/296.raw fake-data/297.raw fake-data/298.raw fake-data/299.raw fake-data/300.raw fake-data/301.raw fake-data/302.raw fake-data/303.raw fake-data/304.raw fake-data/305.raw fake-data/306.raw fake-data/307.raw fake-data/308.raw fake-data/309.raw fake-data/310.raw fake-data/311.raw fake-data/312.raw fake-data/313.raw fake-data/314.raw fake-data/315.raw fake-data/316.raw fake-data/317.raw fake-data/318.raw fake-data/319.raw fake-data/320.raw fake-data/321.raw fake-data/322.raw fake-data/323.raw fake-data/324.raw fake-data/325.raw fake-data/326.raw fake-data/327.raw fake-data/328.raw fake-data/329.raw fake-data/330.raw fake-data/331.raw fake-data/332.raw fake-data/333.raw fake-data/334.raw fake-data/335.raw fake-data/336.raw fake-data/337.raw fake-data/338.raw fake-data/339.raw fake-data/340.raw fake-data/341.raw fake-data/342.raw fake-data/343.raw fake-data/344.raw fake-data/345.raw fake-data/346.raw fake-data/347.raw fake-data/348.raw fake-data/349.raw fake-data/350.raw fake-data/351.raw fake-data/352.raw fake-data/353.raw fake-data/354.raw fake-data/355.raw fake-data/356.raw fake-data/357.raw fake-data/358.raw fake-data/359.raw fake-data/360.raw fake-data/361.raw fake-data/362.raw fake-data/363.raw fake-data/364.raw fake-data/365.raw fake-data/366.raw fake-data/367.raw fake-data/368.raw fake-data/369.raw fake-data/370.raw fake-data/371.raw fake-data/372.raw fake-data/373.raw fake-data/374.raw fake-data/375.raw fake-data/376.raw fake-data/377.raw fake-data/378.raw fake-data/379.raw fake-data/380.raw fake-data/381.raw fake-data/382.raw fake-data/383.raw fake-data/384.raw fake-data/385.raw fake-data/386.raw fake-data/387.raw fake-data/388.raw fake-data/389.raw fake-data/390.raw fake-data/391.raw fake-data/392.raw fake-data/393.raw fake-data/394.raw fake-data/395.raw fake-data/396.raw fake-data/397.raw fake-data/398.raw fake-data/399.raw fake-data/400.raw fake-data/401.raw fake-data/402.raw fake-data/403.raw fake-data/404.raw fake-data/405.raw fake-data/406.raw fake-data/407.raw fake-data/408.raw fake-data/409.raw fake-data/410.raw fake-data/411.raw fake-data/412.raw fake-data/413.raw fake-data/414.raw fake-data/415.raw fake-data/416.raw fake-data/417.raw fake-data/418.raw fake-data/419.raw fake-data/420.raw fake-data/421.raw fake-data/422.raw fake-data/423.raw fake-data/424.raw fake-data/425.raw fake-data/426.raw fake-data/427.raw fake-data/428.raw fake-data/429.raw fake-data/430.raw fake-data/431.raw fake-data/432.raw fake-data/433.raw fake-data/434.raw fake-data/435.raw fake-data/436.raw fake-data/437.raw fake-data/438.raw fake-data/439.raw fake-data/440.raw fake-data/441.raw fake-data/442.raw fake-data/443.raw fake-data/444.raw fake-data/445.raw fake-data/446.raw fake-data/447.raw fake-data/448.raw fake-data/449.raw fake-data/450.raw fake-data/451.raw fake-data/452.raw fake-data/453.raw fake-data/454.raw fake-data/455.raw fake-data/456.raw fake-data/457.raw fake-data/458.raw fake-data/459.raw fake-data/460.raw fake-data/461.raw fake-data/462.raw fake-data/463.raw fake-data/464.raw fake-data/465.raw fake-data/466.raw fake-data/467.raw fake-data/468.raw fake-data/469.raw fake-data/470.raw fake-data/471.raw fake-data/472.raw fake-data/473.raw fake-data/474.raw fake-data/475.raw fake-data/476.raw fake-data/477.raw fake-data/478.raw fake-data/479.raw fake-data/480.raw fake-data/481.raw fake-data/482.raw fake-data/483.raw fake-data/484.raw fake-data/485.raw fake-data/486.raw fake-data/487.raw fake-data/488.raw fake-data/489.raw fake-data/490.raw fake-data/491.raw fake-data/492.raw fake-data/493.raw fake-data/494.raw fake-data/495.raw fake-data/496.raw fake-data/497.raw fake-data/498.raw fake-data/499.raw fake-data/500.raw fake-data/501.raw fake-data/502.raw fake-data/503.raw fake-data/504.raw fake-data/505.raw fake-data/506.raw fake-data/507.raw fake-data/508.raw fake-data/509.raw fake-data/510.raw fake-data/511.raw fake-data/512.raw fake-data/513.raw fake-data/514.raw fake-data/515.raw fake-data/516.raw fake-data/517.raw fake-data/518.raw fake-data/519.raw fake-data/520.raw fake-data/521.raw fake-data/522.raw fake-data/523.raw fake-data/524.raw fake-data/525.raw fake-data/526.raw fake-data/527.raw fake-data/528.raw fake-data/529.raw fake-data/530.raw fake-data/531.raw fake-data/532.raw fake-data/533.raw fake-data/534.raw fake-data/535.raw fake-data/536.raw fake-data/537.raw fake-data/538.raw fake-data/539.raw fake-data/540.raw fake-data/541.raw fake-data/542.raw fake-data/543.raw fake-data/544.raw fake-data/545.raw fake-data/546.raw fake-data/547.raw fake-data/548.raw fake-data/549.raw fake-data/550.raw fake-data/551.raw fake-data/552.raw fake-data/553.raw fake-data/554.raw fake-data/555.raw fake-data/556.raw fake-data/557.raw fake-data/558.raw fake-data/559.raw fake-data/560.raw fake-data/561.raw fake-data/562.raw fake-data/563.raw fake-data/564.raw fake-data/565.raw fake-data/566.raw fake-data/567.raw fake-data/568.raw fake-data/569.raw fake-data/570.raw fake-data/571.raw fake-data/572.raw fake-data/573.raw fake-data/574.raw fake-data/575.raw fake-data/576.raw fake-data/577.raw fake-data/578.raw fake-data/579.raw fake-data/580.raw fake-data/581.raw fake-data/582.raw fake-data/583.raw fake-data/584.raw fake-data/585.raw fake-data/586.raw fake-data/587.raw fake-data/588.raw fake-data/589.raw fake-data/590.raw fake-data/591.raw fake-data/592.raw fake-data/593.raw fake-data/594.raw fake-data/595.raw fake-data/596.raw fake-data/597.raw fake-data/598.raw fake-data/599.raw fake-data/600.raw fake-data/601.raw fake-data/602.raw fake-data/603.raw fake-data/604.raw fake-data/605.raw fake-data/606.raw fake-data/607.raw fake-data/608.raw fake-data/609.raw fake-data/610.raw fake-data/611.raw fake-data/612.raw fake-data/613.raw fake-data/614.raw fake-data/615.raw fake-data/616.raw fake-data/617.raw fake-data/618.raw fake-data/619.raw fake-data/620.raw fake-data/621.raw fake-data/622.raw fake-data/623.raw fake-data/624.raw fake-data/625.raw fake-data/626.raw fake-data/627.raw fake-data/628.raw fake-data/629.raw fake-data/630.raw fake-data/631.raw fake-data/632.raw fake-data/633.raw fake-data/634.raw fake-data/635.raw fake-data/636.raw fake-data/637.raw fake-data/638.raw fake-data/639.raw fake-data/640.raw fake-data/641.raw fake-data/642.raw fake-data/643.raw fake-data/644.raw fake-data/645.raw fake-data/646.raw fake-data/647.raw fake-data/648.raw fake-data/649.raw fake-data/650.raw fake-data/651.raw fake-data/652.raw fake-data/653.raw fake-data/654.raw fake-data/655.raw fake-data/656.raw fake-data/657.raw fake-data/658.raw fake-data/659.raw fake-data/660.raw fake-data/661.raw fake-data/662.raw fake-data/663.raw fake-data/664.raw fake-data/665.raw fake-data/666.raw fake-data/667.raw fake-data/668.raw fake-data/669.raw fake-data/670.raw fake-data/671.raw fake-data/672.raw fake-data/673.raw fake-data/674.raw fake-data/675.raw fake-data/676.raw fake-data/677.raw fake-data/678.raw fake-data/679.raw fake-data/680.raw fake-data/681.raw fake-data/682.raw fake-data/683.raw fake-data/684.raw fake-data/685.raw fake-data/686.raw fake-data/687.raw fake-data/688.raw fake-data/689.raw fake-data/690.raw fake-data/691.raw fake-data/692.raw fake-data/693.raw fake-data/694.raw fake-data/695.raw fake-data/696.raw fake-data/697.raw fake-data/698.raw fake-data/699.raw fake-data/700.raw fake-data/701.raw fake-data/702.raw fake-data/703.raw fake-data/704.raw fake-data/705.raw fake-data/706.raw fake-data/707.raw fake-data/708.raw fake-data/709.raw fake-data/710.raw fake-data/711.raw fake-data/712.raw fake-data/713.raw fake-data/714.raw fake-data/715.raw fake-data/716.raw fake-data/717.raw fake-data/718.raw fake-data/719.raw fake-data/720.raw fake-data/721.raw fake-data/722.raw fake-data/723.raw fake-data/724.raw fake-data/725.raw fake-data/726.raw fake-data/727.raw fake-data/728.raw fake-data/729.raw fake-data/730.raw fake-data/731.raw fake-data/732.raw fake-data/733.raw fake-data/734.raw fake-data/735.raw fake-data/736.raw fake-data/737.raw fake-data/738.raw fake-data/739.raw fake-data/740.raw fake-data/741.raw fake-data/742.raw fake-data/743.raw fake-data/744.raw fake-data/745.raw fake-data/746.raw fake-data/747.raw fake-data/748.raw fake-data/749.raw fake-data/750.raw fake-data/751.raw fake-data/752.raw fake-data/753.raw fake-data/754.raw fake-data/755.raw fake-data/756.raw fake-data/757.raw fake-data/758.raw fake-data/759.raw fake-data/760.raw fake-data/761.raw fake-data/762.raw fake-data/763.raw fake-data/764.raw fake-data/765.raw fake-data/766.raw fake-data/767.raw fake-data/768.raw fake-data/769.raw fake-data/770.raw fake-data/771.raw fake-data/772.raw fake-data/773.raw fake-data/774.raw fake-data/775.raw fake-data/776.raw fake-data/777.raw fake-data/778.raw fake-data/779.raw fake-data/780.raw fake-data/781.raw fake-data/782.raw fake-data/783.raw fake-data/784.raw fake-data/785.raw fake-data/786.raw fake-data/787.raw fake-data/788.raw fake-data/789.raw fake-data/790.raw fake-data/791.raw fake-data/792.raw fake-data/793.raw fake-data/794.raw fake-data/795.raw fake-data/796.raw fake-data/797.raw fake-data/798.raw fake-data/799.raw fake-data/800.raw fake-data/801.raw fake-data/802.raw fake-data/803.raw fake-data/804.raw fake-data/805.raw fake-data/806.raw fake-data/807.raw fake-data/808.raw fake-data/809.raw fake-data/810.raw fake-data/811.raw fake-data/812.raw fake-data/813.raw fake-data/814.raw fake-data/815.raw fake-data/816.raw fake-data/817.raw fake-data/818.raw fake-data/819.raw fake-data/820.raw fake-data/821.raw fake-data/822.raw fake-data/823.raw fake-data/824.raw fake-data/825.raw fake-data/826.raw fake-data/827.raw fake-data/828.raw fake-data/829.raw fake-data/830.raw fake-data/831.raw fake-data/832.raw fake-data/833.raw fake-data/834.raw fake-data/835.raw fake-data/836.raw fake-data/837.raw fake-data/838.raw fake-data/839.raw fake-data/840.raw fake-data/841.raw fake-data/842.raw fake-data/843.raw fake-data/844.raw fake-data/845.raw fake-data/846.raw fake-data/847.raw fake-data/848.raw fake-data/849.raw fake-data/850.raw fake-data/851.raw fake-data/852.raw fake-data/853.raw fake-data/854.raw fake-data/855.raw fake-data/856.raw fake-data/857.raw fake-data/858.raw fake-data/859.raw fake-data/860.raw fake-data/861.raw fake-data/862.raw fake-data/863.raw fake-data/864.raw fake-data/865.raw fake-data/866.raw fake-data/867.raw fake-data/868.raw fake-data/869.raw fake-data/870.raw fake-data/871.raw fake-data/872.raw fake-data/873.raw fake-data/874.raw fake-data/875.raw fake-data/876.raw fake-data/877.raw fake-data/878.raw fake-data/879.raw fake-data/880.raw fake-data/881.raw fake-data/882.raw fake-data/883.raw fake-data/884.raw fake-data/885.raw fake-data/886.raw fake-data/887.raw fake-data/888.raw fake-data/889.raw fake-data/890.raw fake-data/891.raw fake-data/892.raw fake-data/893.raw fake-data/894.raw fake-data/895.raw fake-data/896.raw fake-data/897.raw fake-data/898.raw fake-data/899.raw fake-data/900.raw fake-data/901.raw fake-data/902.raw fake-data/903.raw fake-data/904.raw fake-data/905.raw fake-data/906.raw fake-data/907.raw fake-data/908.raw fake-data/909.raw fake-data/910.raw fake-data/911.raw fake-data/912.raw fake-data/913.raw fake-data/914.raw fake-data/915.raw fake-data/916.raw fake-data/917.raw fake-data/918.raw fake-data/919.raw fake-data/920.raw fake-data/921.raw fake-data/922.raw fake-data/923.raw fake-data/924.raw fake-data/925.raw fake-data/926.raw fake-data/927.raw fake-data/928.raw fake-data/929.raw fake-data/930.raw fake-data/931.raw fake-data/932.raw fake-data/933.raw fake-data/934.raw fake-data/935.raw fake-data/936.raw fake-data/937.raw fake-data/938.raw fake-data/939.raw fake-data/940.raw fake-data/941.raw fake-data/942.raw fake-data/943.raw fake-data/944.raw fake-data/945.raw fake-data/946.raw fake-data/947.raw fake-data/948.raw fake-data/949.raw fake-data/950.raw fake-data/951.raw fake-data/952.raw fake-data/953.raw fake-data/954.raw fake-data/955.raw fake-data/956.raw fake-data/957.raw fake-data/958.raw fake-data/959.raw fake-data/960.raw fake-data/961.raw fake-data/962.raw fake-data/963.raw fake-data/964.raw fake-data/965.raw fake-data/966.raw fake-data/967.raw fake-data/968.raw fake-data/969.raw fake-data/970.raw fake-data/971.raw fake-data/972.raw fake-data/973.raw fake-data/974.raw fake-data/975.raw fake-data/976.raw fake-data/977.raw fake-data/978.raw fake-data/979.raw fake-data/980.raw fake-data/981.raw fake-data/982.raw fake-data/983.raw fake-data/984.raw fake-data/985.raw fake-data/986.raw fake-data/987.raw fake-data/988.raw fake-data/989.raw fake-data/990.raw fake-data/991.raw fake-data/992.raw fake-data/993.raw fake-data/994.raw fake-data/995.raw fake-data/996.raw fake-data/997.raw fake-data/998.raw fake-data/999.raw

grind: unit_library
	valgrind --leak-check=full --error-exitcode=2 ./unit_library -quiet

clean:
	rm -f unit_geometry unit_lists unit_storage unit_image unit_align unit_blob unit_model unit_library
	rm -f test_image.tiff performance_imprint.tiff worm_imprint.tiff worm_noisy.tiff
	rm -f test_blob.log
	rm -f *.o
	rm -f test*.tiff
	rm -f 20071226_105033/*
	rm -f 20071212_130514/*
	if [ -a 20071226_105033 ]; then rmdir 20071226_105033; fi
	if [ -a 20071212_130514 ]; then rmdir 20071212_130514; fi
	rm -f lib/MWT.dll
	if [ -a lib ]; then rmdir lib; fi
