import subprocess
from collections import Counter

error_lines = """53289699            800
53289724           4035
53289732           5248
53289762            373
53289947           4035
53289959           5507
53289961           5513
53289985           7857
53289994           9314
53290060          17338
53290096          19928
53290165          25152
53290186          25368
53290201          27059
53290224          30975
53290254          36420
53290255          36422
53290272          41535
53290278          41910
53290279          42080
53290281          42088
53290284          43366
53290285          43370
53290286          43434
53290287          43827
53290298          47429
53290299          47431
53290300          47433
53290343          51930
53290344          51939
53290367          55074
53290375          58445
53290376          58704
53290380          59630
53290391          61131
53290406          63168
53290418          64960
53290451          71729
53290452          71853
53290610          73717
53290833          76664
53290834          76665
53291290          85096
53291294          85261
53291296          85568
53291298          85571
53291308          87502
53291310          88394
53291312          88704
53291323          92477
53291333          95508
53291340          96925
53291342          96969
53291354          97885
53291359          98065
53291364          99644
53291372         101568
53291390         107549
53291421         111892
53291449         117394
53291450         117396
53291451         117397
53291463         121536
53291467         122580
53291472         122625
53291473         122626
53291474         122627
53291490         127374
53291495         127913
53291511         131334
53291516         132717
53291519         132769
53291522         133410
53291525         135485
53291526         135491
53291534         137847
53291535         137852
53291537         138156
53291538         138157
53291555         139424
53291558         139829
53291565         142882
53291569         143681
53291573         144081
53291577         153156
53291578         153171
53291579         153770
53291582         154528
53291586         155059
53291587         155412
53291588         155413
53291593         155571
53291620         157085
53291673         158641
53291686         159552
53291696         159555
53291727         161173
53291757         163340
53291759         163398
53291764         163402
53291813         168024
53291894         178249
53291899         179089
53291904         180530
53291911         181015
53291912         181017
53291915         181025
53291918         181026
53291920         182210
53291949         189227
53291950         189229
53291951         189232
53291953         190589
53291954         190724
53291988         197714
53291992         197944
53291998         198418
53292005         198793
53292009         199512
53292010         199513
53292011         199515
53292012         199517
53292079         200763
53292081         201113
53292082         201268
53292276         204361
53292277         204365
53292278         204368
53292341         205057
53292538         206742
53292602         208433
53292603         208435
53292737         209257
53292870         211579
53293067         213747
53293075         214168
53293114         219160
53293120         221534
53293134         222484
53293141         224548
53293142         224549
53293153         225011
53293178         228721
53293188         229475
53293194         229952
53293204         231406
53293207         231496
53293225         235393
53293237         237050
53293243         238712
53293249         239321
53293251         239666
53293252         239667
53293253         239673
53293261         242442
53293262         242454
53293269         243562
53293272         243576
53293280         243953
53293286         244089
53293287         244090
53293302         245202
53293303         245208
53293304         245209
53293305         245216
53293306         245217
53293307         245218
53293308         245224
53293309         245225
53293310         245232
53293311         245234
53293312         245235
53293315         245241
53293325         245266
53293326         245272
53293328         245280
53293330         245287
53293331         245289
53293332         245291
53293333         245298
53293334         245299
53293335         245304
53293336         245305
53293337         245312
53293339         245313
53293340         245322
53293341         245323
53293342         245328
53293343         245330
53293344         245331
53293345         245333
53293348         245336
53293349         245337
53293350         245340
53293352         245344
53293353         245345
53293355         245352
53293356         245353
53293357         245360
53293358         245361
53293360         245368
53293361         245371
53293364         245380
53293365         245384
53293368         245412
53293369         245414
53293370         245416
53293378         245947
53293384         245991
53293388         245996
53293390         246002
53293391         246017
53293398         246219
53293411         247195
53293413         247197
53293434         251355
53293445         252367
53293448         253418
53293499         253892
53293519         255902
53293521         255907
53293529         256338
53293531         257042
53293534         257256
53293536         257258
53293545         258552
53293561         261524
53293566         262406
53293575         262659
53293578         262690
53293590         263779
53293624         267749
53293665         269073
53293669         269075
53293672         269077
53293676         269079
53293720         271398
53293752         272242
53293788         274816
53293808         276651
53293811         276941
53293829         277778
53294210            800
53294247           3371
53294254           4035
53294295           5254
53294347           7857
53294366           9314
53294372           9601
53294413          15376
53294414          15377
53294431          17211
53294472          19928
53294481          20169
53294510          24694
53294598          25368
53294603          27059
53294756          30975
53294850          33970
53295110          35408
53295126          36420
53295128          36422
53295165          41187
53295169          41535
53295171          41723
53295176          41910
53295177          42080
53295178          42088
53295181          43366
53295182          43370
53295183          43434
53295187          43847
53295222          47429
53295226          47431
53295229          47433
53295293          51378
53295383          55074
53295391          56101
53295409          58445
53295413          58704
53295541          64960
53295559          68126
53295574          71336
53295576          71580
53295590          72545
53295607          73717
53295623          76664
53295625          76665
53295675          85096
53295717          88374
53295726          88704
53295732          88819
53295757          92477
53295780          95508
53295795          96969
53295807          97885
53295820          99644
53295822          99646
53295829         101568
53295840         103489
53295851         107549
53295944         111892
53296025         122625
53296027         122627
53296033         123213
53296043         126204
53296051         127374
53296097         131334
53296112         132717
53296119         132769
53296124         133410
53296127         135485
53296128         135491
53296136         137847
53296137         137852
53296139         138156
53296140         138157
53296161         139424
53296163         139829
53296168         142882
53296173         143681
53296199         144081
53296202         153156
53296203         153171
53296204         153770
53296207         154528
53296211         155058
53296214         155059
53296216         155412
53296219         155413
53296227         155571
53296260         157085
53296287         158213
53296291         158216
53296300         158641
53296304         159552
53296307         159555
53296324         161173
53296336         163340
53296340         163398
53296343         163402
53296379         170102
53296395         171511
53296417         176422
53296457         179089
53296470         181015
53296471         181017
53296474         181025
53296479         182210
53296540         187923
53296560         190589
53296561         190724
53296590         197376
53296592         197714
53296594         197944
53296598         198422
53296603         198793
53296606         199512
53296608         199513
53296609         199515
53296610         199517
53296622         200763
53296630         201268
53296679         204361
53296684         204365
53296686         204368
53296699         205057
53296728         209257
53296819         213747
53296826         214168
53296845         221534
53296854         222293
53296857         222484
53296904         228721
53296917         229475
53296923         229952
53296933         231363
53296938         231496
53296945         233476
53296975         237050
53296984         238712
53296987         239150
53296990         239321
53296994         239667
53296995         239673
53297003         242442
53297004         242454
53297012         243562
53297015         243576
53297031         244084
53297032         244089
53297048         245202
53297050         245209
53297054         245224
53297099         245337
53297102         245345
53297105         245353
53297107         245361
53297120         245418
53297128         245952
53297198         251355
53297236         253892
53297269         257042
53297271         257256
53297273         257258
53297277         258518
53297279         258532
53297280         258549
53297281         258552
53297282         258553
53297299         262406
53297308         262659
53297311         262690
53297337         267748
53297338         267749
53297350         269073
53297352         269075
53297353         269077
53297354         269079
53297363         270476
53297375         272242
53297385         274816
53297391         276543
53297393         276651
53297394         276941
53297400         277778
53297414         278191
53297447         282229
53297448         282231
53297449         282235
53297469         285265
53297471         285270
53297472         285271
53297473         285474
53297480         286504
53297488         288301
53297489         288303
53297490         288304
53297491         288306
53297492         288308
53297495         289257
53297516         291710
53297531         292880
53297534         292889
53297536         293843
53297537         293845
53297550         295984
53297569         299032
53297573         300546
53297576         301193
53297578         301211
53297579         301217
53297676         303315
53297681         303317
53297685         303357
53297749         307659
53297751         308057
53297753         308058
53297772         308120
53297773         308170
53297808         310685
53297821         312167
53297824         313057
53297846         315818
53297847         315833
53297854         316929
53297855         317430
53297870         321026
53297872         321028
53297883         326059
53297886         328384
53297934         335343
53297990         338608
53298042         343554
53298082         345033
53298187         353660
53298209         359154
53298251         362331
53298270         364380
53298272         364382
53298285         364843
53298286         364849
53298287         364850
53298288         364851
53298301         366030
53298322         368900
53298330         371798
53298331         371799
53298333         371801
53298335         371949
53298339         372083
53298349         374329
53298353         374885
53298354         375276
53298404         378564
53298446         380467
53298451         380468
53298453         380472
53298487         381911
53298490         381916
53298496         381917
53298501         381918
53298504         381921
53298506         381927
53298518         383280
53298521         383281
53298523         383284
53298525         383288
53298534         383312
53298539         383322
53298543         383328
53298548         383329
53298596         387222
53298627         389385
53298644         393932
53298646         393936
53298648         393937
53298670         395424
53298671         395426
53298674         395861
53298675         395864
53298676         395865
53298677         395866
53298678         395891
53298680         395894
53298684         397366
53298685         397373
53298688         398485
53298689         398502
53298690         398592
53298691         398593
53298692         398597
53298734         401532
53298737         401547
53298738         401694
53298739         401695
53298740         401696
53298741         401697
53298742         401962
53298750         402206
53298751         402705
53298775         407542
53298780         407845
53298781         407846
53298792         411858
53298811         414709
53298812         414830
53298814         414833
53298815         415708
53298820         416121
53298856         420595
53298858         420599
53298865         421142
53298867         421178
53298907         430630
53298908         430631
53298921         432232
53298969         439663
53298970         439683
53298971         439686
53298975         441193
53298978         441283
53298979         441422
53298988         442605
53298989         442611
53298996         442765
53299000         442771
53299001         442776
53299011         443979
53299031         446760
53299033         446785
53299085         453466
53299092         454569
53299095         454978
53299108         458305
53299109         458307
53299120         458915
53299132         460677
53299138         461363
53299170         463179
53299177         463727
53299197         466253
53299219         470343
53299220         470408
53299221         470416
53304192           3361
53304232           4035
53304480           7857
53304536           9314
53304603          16309
53304615          18888
53304625          19928
53304629          20159
53304656          24694
53304671          25368
53304682          27059
53304703          30929
53304704          30975
53304753          36420
53304754          36422
53304770          41535
53304777          41910
53304778          42080
53304780          42088
53304783          43366
53304786          43370
53304787          43434
53304800          47429
53304803          47431
53304804          47433
53304877          51368
53304893          52467
53304908          55074
53304916          58445
53304917          58704
53304923          59637
53304971          64955
53304972          64960
53304974          64961
53305064          71711
53305078          72545
53305080          72547
53305089          73717
53305101          75675
53305132          76664
53305134          76665
53305195          85567
53305213          88704
53305227          92477
53305237          95508
53305255          97360
53305261          97885
53305274          99644
53305284         101568
53305301         107549
53305361         111892
53305384         116929
53305417         121538
53305435         122625
53305437         122627
53305442         123213
53305447         126197
53305452         127374
53305458         127927
53305472         131334
53305474         131653
53305477         132717
53305480         132769
53305483         133410
53305486         135485
53305487         135491
53305494         137846
53305495         137847
53305496         137852
53305498         138156
53305499         138157
53305519         139424
53305521         139829
53305526         142882
53305530         143681
53305533         144081
53305536         153156
53305537         153171
53305538         153770
53305540         154528
53305545         155059
53305546         155412
53305547         155413
53305551         155571
53305561         157085
53305577         158641
53305582         159552
53305585         159555
53305595         161173
53305599         161297
53305622         163340
53305623         163398
53305624         163402
53305675         176422
53305687         178248
53305693         179089
53305706         181015
53305719         181025
53305720         181026
53305726         182210
53305732         183097
53305737         183554
53305765         190589
53305766         190724
53305794         197714
53305796         197944
53305809         198793
53305812         199512
53305813         199513
53305814         199515
53305815         199517
53305816         199536
53305821         200763
53305824         201268
53305838         204361
53305839         204365
53305840         204368
53305844         205057
53305860         208649
53305865         209257
53305885         213747
53305909         219160
53305913         221534
53305925         222484
53305972         228721
53305987         229475
53305995         229952
53306010         231496
53306041         237050
53306047         238712
53306053         239321
53306056         239667
53306057         239673
53306065         242442
53306066         242454
53306073         243562
53306076         243576
53306125         245249
53306139         245289
53306140         245291
53306154         245337
53306238         251355
53306288         253892
53306310         255902
53306323         257042
53306325         257256
53306327         257258
53306333         258532
53306334         258549
53306335         258552
53306337         258553
53306353         260166
53306361         262406
53306369         262659
53306372         262690
53306399         267749
53306436         269073
53306437         269075
53306438         269077
53306439         269079
53306451         270480
53306461         271626
53306476         272242
53306488         274816
53306494         276543
53306496         276651
53306497         276941
53306503         277778
53306560         282229
53306561         282231
53306562         282235
53306582         285265
53306584         285270
53306585         285271
53306586         285474
53306602         288301
53306603         288303
53306604         288304
53306605         288306
53306607         288308
53306610         289257
53306634         291710
53306672         292880
53306676         292889
53306678         293843
53306679         293845
53306694         295984
53306715         299032
53306718         300546
53306720         301193
53306721         301197
53306722         301211
53306723         301217
53306726         301223
53306732         301230
53306737         301508
53306738         301511
53306756         303357
53306774         306522
53306832         307659
53306834         308057
53306835         308058
53306836         308074
53306837         308076
53306838         308090
53306839         308092
53306841         308120
53306842         308170
53306858         308960
53306892         310685
53306894         310808
53306909         313057
53306917         313737
53306929         315818
53306930         315833
53306937         316929
53306947         319216
53306948         319222
53306954         321026
53306958         321028
53306988         326059
53306994         328384
53307020         333902
53307044         336014
53307066         338608
53307071         338868
53307083         341984
53307116         345033
53307117         345082
53307120         345129
53307121         345130
53307122         345132
53307176         353660
53307216         359154
53307289         362331
53307350         366030
53307359         367890
53307371         368900
53307378         371781
53307385         371794
53307387         371799
53307388         371800
53307393         371949
53307398         372083
53307410         374329
53307413         374885
53307417         375276
53307498         378564
53307555         380467
53307561         380468
53307566         380472
53307635         381911
53307640         381916
53307645         381917
53307651         381918
53307657         381921
53307663         381927
53307687         383280
53307692         383281
53307696         383284
53307702         383288
53307711         383312
53307716         383322
53307722         383328
53307728         383329
53307774         387222
53307794         389385
53307827         393932
53307828         393936
53307829         393937
53307838         394369
53307839         394455
53307846         395424
53307847         395426
53307851         395861
53307852         395864
53307853         395865
53307856         395866
53307857         395891
53307858         395894
53307860         397366
53307861         397373
53307864         398485
53307866         398502
53307867         398592
53307868         398593
53307869         398597
53307916         401547
53307917         401694
53307918         401695
53307919         401696
53307920         401697
53307922         401962
53307929         402705
53307930         402707
53307931         402719
53307932         402723
53307950         407070
53307954         407542
53307956         407664
53307969         411858
53307987         414417
53307989         414709
53307990         414830
53307991         414833
53308026         419806
53308045         421142
53308046         421178
53308074         427841
53308088         430631
53308164         439663
53308165         439683
53308166         439686
53308171         441193
53308173         441283
53308174         441422
53308175         441423
53308188         442765
53308190         442771
53308199         443979
53308206         444478
53308218         446760
53308219         446785
53308263         453466
53308289         454569
53308312         458872
53308315         458915
53308323         460677
53308327         461363
53308334         462130
53308341         463179
53308348         463727
53308363         466269
53308373         469289
53308374         469505
53308382         470343
53308383         470408
53308384         470416
53308409         475500
53308440         479050
53308448         479681
53308456         482285
53308948         491265
53308950         491266
53308951         491267
53308952         491268
53308953         491269
53308954         491270
53308955         491400
53308980         493512
53308989         494069
53309015         497883
53309028         500407
53309309         510383
53309358         511550
53309426         513034
53309434         513311
53309521         515004
53309523         515165
53309562         516554
53309565         517008
53309588         517877
53309593         518300
53309596         518304
53309597         518306
53309598         518308
53309600         518312
53309601         518314
53309609         518576
53309616         518578
53309622         518588
53309629         518592
53309675         520145
53309676         520149
53309687         521095
53309688         521155
53309689         521188
53309690         521221
53309691         521227
53309695         521228
53309696         521255
53309706         522608
53309792         525348
53309793         525372
53309811         526490
53309812         526491
53309833         528152
53309834         528153
53309841         529164
53309850         529974
53309858         530193
53309865         530725
53309866         530735
53309867         530736
53309868         530737
53309869         530739
53310076         532266
53310085         533206
53310091         533690
53310092         533691
53310093         533700
53310094         533701
53310100         534704
53310111         535236
53310121         536317
53310122         536323
53310124         536625
53310125         536768
53310126         536769
53310127         536832
53310131         537190
53310137         537545
53310153         538664
53310155         538677
53310156         538678
53310157         538679
53310159         538706
53310161         538718
53310162         538719
53310163         538793
53310166         539109
53310167         539179
53310193         539836
53310198         540425
53310200         540436
53310201         540437
53310202         540438
53310216         540460
53310217         540464
53310219         540466
53310221         540468
53310223         540470
53310262         543710
53310263         544110
53310264         544112
53310265         544116
53310266         544118
53310288         545418
53310293         545420
53310307         545530
53310308         545545
53310315         546304
53310317         546461
53310318         546462
53310324         547221
53310325         547329
53310326         547330
53310327         547337
53310328         547341
53310329         547342
53310330         547348
53310331         547361
53310332         547392
53310336         547706
53310337         547707
53310351         548009
53310354         548011
53310390         549528
53310399         550848
53310401         550912
53310419         554361
53310422         554791
53310424         554931
53310425         554937
53310427         554939
53310428         554941
53310429         555200
53310737         562457
53310739         562463
53310755         563506
53310757         563763
53310759         563768
53310862         563770
53311246         564614
53311283         566088
53311301         566096
53311331         567263
53311500         572417
53311529         573842
53311531         573846
53311574         575814
53311601         577491
53311602         577492
53311603         577493
53311604         577494
53311605         577495
53311606         577496
53311607         577499
53311608         577500
53311609         577501
53311610         577502
53311611         577503
53311627         580088
53311667         585023
53311692         588564
53311695         588567
53311732         594094
53311736         595521
53311760         598164
53311786         599869
53311788         600194
53311789         600196
53311801         601397
53311802         601399
53311805         602139
53311806         602140
53311807         602142
53311816         603706
53311822         604593
53311832         608735
53311836         609175
53311864         613941
53311866         613954
53311867         613955
53311869         613977
53311870         613993
53311871         613999
53311881         616051
53311889         616871
53311898         618188
53311935         622518
53311945         623544
53311961         624427
53311985         625642
53312013         628903
53312045         630536
53312056         630962
53312065         632020
53312073         632519
53312087         633766
53312088         633767
53312089         633768
53312099         635058
53312119         636915
53312120         637002
53312121         637003
53312124         637256
53312129         637315
53312139         638414
53312226         640887
53312237         641374
53312262         642737
53312287         643021
53312288         643025
53312289         643027
53312290         643029
53312291         643033
53312292         643035
53312342         649288
53312347         649790
53312361         651506
53312363         651508
53312372         653335
53312383         654985
53312389         656385
53312391         656438
53312398         656951
53312399         656975
53312400         656976
53312405         660144
53312416         661888
53312417         661890
53312418         661891
53312421         662018
53312422         662272
53312423         662786
53312424         662787
53312425         663852
53312426         663882
53312427         663936
53312428         664203
53312429         664513
53312430         664705
53312431         665346
53312432         665632
53312433         665729
53312439         665981
53312440         665991
53312465         667841
53312466         667873
53312467         667885
53312468         667888
53312469         667911
53312470         667968
53312471         668107
53312472         668111
53312473         668156
53312474         668192
53312475         668193
53312476         668320
53312479         668691
53312484         668932
53312486         668934
53312500         671413"""
ids = subprocess.check_output("sacct -S 11/07-17:10 | grep FAILED | grep norm | cut -c-24", shell=True).split("\n")
errors = Counter()
for id in ids:
    if not id: continue
    slurm_id, ibis_id = id.split()
    try:
        print "slurm-{}.out".format(slurm_id)
        with open("slurm-{}.out".format(slurm_id)) as f:
            for last_line in f: pass
            error = last_line.rstrip().split(":")[0]
            errors[error] += 1

            if error  in [
              #"struct.error",
              #"KeyError",
              #"cElementTree.ParseError",
              #"OSError",
              #"biopdbtools.InvalidPDB",
              "IOError"
              ]:
                with open("{}.sh".format(ibis_id)) as f2:
                    for last_line_i in f2: pass
                print last_line_i
                f.seek(0)
                print f.read()
    except IOError as e:
        print e
print
for error, count in errors.iteritems():
    print error, count
