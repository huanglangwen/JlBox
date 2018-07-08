const RO2_names=["HOCH2CH2O2","HO1C3O2","HYPROPO2","IPROPOLO2","NBUTOLAO2","NBUTOLBO2","BUT2OLO2","IBUTOLBO2","IBUTOLCO2","TBUTOLO2","HO3C5O2",
"PE2ENEBO2","HM2C43O2","M2BUOL2O2","HM33C3O2","ME3BUOLO2","HO2M2C4O2","ME2BU2OLO2","PROL11MO2","H2M3C4O2","ME2BUOLO2","CYHEXOLAO2",
"MIBKAOHAO2","MIBKAOHBO2","MIBKHO4O2","CH3CO3","NMBOAO2","NMBOBO2","MBOAO2","MBOBO2","CH3O2","HCOCH2O2","C2H5CO3","C2H5O2","C3H7CO3",
"BUTALO2","NC3H7O2","IPRCO3","IBUTALBO2,""IBUTALCO2","IC3H7O2","C4H9CO3","C4CHOBO2","NC4H9O2","ACO3","ACRO2","OCCOHCO2","CH3C2H2O2",
"MACO3","MACRO2","MACROHO2","C3DBCO3","C4CONO3O2","C4NO3COO2","C4OCCOHCO2","COCCOH2CO2","SC4H9O2","IC4H9O2","TC4H9O2","PEAO2","PEBO2",
"PECO2","IPEAO2","IPEBO2","IPECO2","NEOPO2","HEXAO2","HEXBO2","HEXCO2","M2PEAO2","M2PEBO2","M2PECO2","M2PEDO2","M3PEAO2","M3PEBO2",
"M3PECO2","M22C43O2","M22C4O2","M33C4O2","M23C43O2","M23C4O2","HEPTO2","M2HEXAO2","M2HEXBO2","M3HEXAO2","M3HEXBO2","OCTO2","NONO2",
"DECO2","UDECO2","DDECO2","CHEXO2","ETHENO3O2","PRONO3AO2","PRONO3BO2","BU1ENO3O2","C43NO34O2","HO3C4O2","C42NO33O2","MPRANO3O2",
"MPRBNO3O2","C51NO32O2","C52NO31O2","PE1ENEAO2","PE1ENEBO2","C52NO33O2","C53NO32O2","PE2ENEAO2","C4NO32M1O2","C4NO32M2O2","C4M3NO31O2",
"C4M3NO32O2","IPRCHOOA","ME3BU2OLO2","C4M2NO32O2","C4M2NO33O2","C65NO36O2","C66NO35O2","C6OH5O2","HO5C6O2","C62NO33O2)","C63NO32O2",
"C64OH5O2","C65OH4O2","C4ME2NO3O2","C4ME2OHO2","BZBIPERO2","C6H5CH2O2","TLBIPERO2","OXYLO2","OXYBIPERO2","MXYLO2","MXYBIPERO2",
"PXYLO2","PXYBIPERO2","C6H5C2H4O2","EBZBIPERO2","PHC3O2","PBZBIPERO2","PHIC3O2","IPBZBIPRO2)","TM123BO2","TM123BPRO2","TM124BO2",
"TM124BPRO2","TMBO2","TM135BPRO2","ETOLO2","OETLBIPRO2","METLBIPRO2","PETLBIPRO2","DM35EBO2","DMEBIPRO2","DE35TO2","DETLBIPRO2",
"NSTYRENO2","STYRENO2","C6H5CO3","C6H5O2","CH2CLO2","CHCL2O2","CCL3O2","CCL3CH2O2","TCEOHO2","C2CL3OHAO2)","C2CL3OHBO2","C2CL2OHO2",
"DICLETO2","CH2OHCL2O2","CL2OHCH2O2","CL12PRAO2","CL12PRBO2","CL12PRCO2","CH3CCL2O2","CHCL2CH2O2","CH2CLCH2O2","CH3CHCLO2",
"CHCL2CL2O2","CH2CL3O2","CHCL3O2","CCLNO3O2)","CNO3CLO2","CCLOHCH2O2","CH2OHCCLO2","NBUTDAO2","NBUTDBO2","BUTDAO2","BUTDBO2",
"BUTDCO2","NISOPO2","ISOP34O2","CHOOCH2O2","METHACETO2","MOCOCH2O2","ACETC2H4O2","EOCOCH2O2","ETHACETO2","NPROACEAO2",
"NPROACEBO2","NPROACECO2","IPRACBO2","IPROACETO2","NBUACETAO2","NBUACETBO2","NBUACETCO2","SBUACETAO2","SBUACETBO2","MCOOTBO2",
"TBOCOCH2O2","CH3OCH2O2","DIETETO2","ETOC2O2","MTBEAO2","MTBEBO2","DIIPRETO2","IPROMC2O2","ETBEAO2","ETBEBO2)","ETBECO2","MO2EOLAO2",
"MO2EOLBO2","EOX2EOLAO2","EOX2EOLBO2","PR2OHMOXO2","H2C3OCO2","BOX2EOHAO2","BOX2EOHBO2","BOXPROLAO2","BOXPROLBO2","CH2BRO2",
"DIBRETO2","CH3COCH2O2","MEKAO2","MEKBO2","MEKCO2","CO2C54O2","MPRKAO2","DIEKAO2","DIEKBO2","MIPKAO2","MIPKBO2","HEX2ONAO2",
"HEX2ONBO2","HEX2ONCO2","HEX3ONAO2","HEX3ONBO2","HEX3ONCO2","HEX3ONDO2","MIBKAO2","MIBKBO2","MTBKO2","CYHXONAO2","NAPINAO2",
"NAPINBO2","APINAO2","APINBO2","APINCO2","NBPINAO2","NBPINBO2","BPINAO2","BPINBO2","BPINCO2","NLIMO2","LIMAO2","LIMBO2","LIMCO2",
"NBCO2","BCAO2","BCBO2","BCCO2","DMMAO2","DMMBO2","DMCO2","CH3SCH2O2","HODMSO2","ETHOXO2","BUT2CO3","C3ME3CO3","C3ME3CHOO2",
"HOCH2CO3)","CH3CHOHCO3","IPRHOCO3","IPRCHOO","BZEMUCCO3","BZEMUCO2","C5DIALO2","NPHENO2","PHENO2","CRESO2","NCRESO2","TLEMUCCO3",
"TLEMUCO2","C615CO2O2","OXYMUCCO3","OXYMUCO2","MC6CO2O2","NOXYOLO2","OXYOLO2","MXYMUCCO3","MXYMUCO2","C726CO5O2","MXYOLO2)","NMXYOLO2",
"PXYMUCCO3","PXYMUCO2","C6M5CO2O2","NPXYOLO2","PXYOLO2","EBENZOLO2","NEBNZOLO2","EBZMUCCO3","EBZMUCO2","C715CO2O2","NPBNZOLO2",
"PBENZOLO2","PBZMUCCO3","PBZMUCO2","C815CO2O2","IPBENZOLO2","NIPBNZOLO2","IPBZMUCCO3","IPGLOOB","IPBZMUCO2","C7M15CO2O2","NTM123OLO2",
"TM123OLO2","TM123MUCO2","NTM124OLO2","TM124OLO2","TM124MUCO3","TM124MUCO2","C7CO2M5O2","NTM135OLO2","TM135OLO2","TM135MUCO3",
"TM135MUCO2","C7M2CO5O2","OETLMUCCO3","OETLMUCO2","MC7CO2O2","NOETOLO2","OETOLO2","METLMUCCO3)","METLMUCO2","C826CO3O2","METOLO2",
"NMETOLO2","PETLMUCCO3","PETLMUCO2","C7M6CO2O2","NPETOLO2","PETOLO2","DMEBMUCO3","DMEBMUCO2","C8M2CO6O2","NDMEPHOLO2","DMEPHOLO2",
"NDEMPHOLO2","DEMPHOLO2","DETLMUCO3","DETLMUCO2","C9M2CO6O2","HMVKAO2","HMVKBO2","MVKO2","CISOPAO2","ISOPBO2","CISOPCO2","ISOPDO2",
"NC526O2","C530O2","M3BU3ECO3","C45O2","NC51O2","C51O2","CH2CHCH2O2","ISOPAO2","ISOPCO2","MEMOXYCO3","EOX2MECO3","ETOMEO2",
"PRONEMOXO2","BOXMCO3","BOX2MO2","BOXPRONAO2","BOXPRONBO2","C107O2","C109O2","C96O2","NOPINAO2","NOPINBO2","NOPINCO2","NOPINDO2",
"LIMALAO2","LIMALBO2","C923O2","BCALAO2","BCALBO2","C136O2","BCALCO2","C141O2","HOC2H4CO3","HOIPRCO3","HO13C5O2","HO3C4CO3","C54O2",
"H2M2C3CO3","PROL1MCO3","C56O2","HO2C43CO3","MIBKCOOHO2","NC4OHCO3","C4OH2CO3","CO2C3CO3","HO2C3CO3","IBUDIALCO3","PROPALO2",
"CO3C4CO3","HO1C4O2","A2PANOO","HCOCOHCO3","HCOCO3","MACRNCO3","MACRNBCO3","CHOMOHCO3","CO2H3CO3","HO1C5O2","HO2C5O2","C52O2",
"TBUTCO3","HO1C6O2","C5H11CO3","HO2C6O2","HO3C6O2","HO1MC5O2","C54CO3","HO2MC5O2","EIPKAO2","EIPKBO2","HO2M2C5O2","H1MC5O2",
"M3C4CO3","H2MC5O2","M2BKAO2","M2BKBO2","HM33C4O2","M22C3CO3","HM22C4O2","M33C3CO3","HM23C4O2","M2C43CO3","HO3C76O2","CO3C75O2",
"H2M5C65O2","C75O2","H2M2C65O2","H2M4C65O2","C710O2","H3M3C6O2","HO3C86O2","CO3C85O2","HO3C96O2","C91O2","HO3C106O2","C101O2",
"HO3C116O2","C111O2","HO3C126O2","C121O2","CO1C6O2","NO3CH2CO3","PRNO3CO3","CO3C4NO3O2","HO3C3CO3","MPRBNO3CO3","C5NO3COAO2",
"C4NO3CO3","C5OH2CO4O2","C4OHCO3","C5NO3CO4O2","C5CONO34O2","C43NO3CO3","C4MCONO3O2","C3MNO3CO3","C3M3OH2CO3","MC4CONO3O2",
"C65NO36CO3","MNO3COC4O2","C4COMOH3O2","HO5C5CO3","C6NO3CO5O2","C6CONO34O2","MALDIALCO3","EPXDLCO3","C3DIALO2","MALDIALO2",
"OXYL1O2","C5CO14O2","OXYLCO3","EPXM2DLCO3","C4MCO2O2","DM123O2","MXYLCO3","MXYL1O2","C3MCODBCO3","EPXMDLCO3","C3MDIALO2",
"MXY1O2","PXYLCO3","PXYL1O2","PXY1O2","C6H5CH2CO3","EBENZO2","C6DCARBBO2","PHCOETO2","PBENZO2","C7DCCO3","IPBENZO2","IC7DCCO3",
"IPGLOO","TM123BCO3","TM123O2","EPXKTMCO3","C4CO2O2","TM124BCO3","DM124O2","TM124O2","TMBCO3","DMPHO2","C4MCODBCO3","EPXMKTCO3",
"CO24C53O2","MPHCOMEO2","EPXMEDLCO3","C4ECO2O2","OET1O2","MET1O2","PET1O2","DMPHCOMO2","EMPHCOMO2","EMPHCO3","C7CODBCO3","EPXEKTCO3",
"C3EDIALO2","CO24C63O2","CCL3CO3","CLETO3","CL2OHCO3","CL12CO3","CLCOCLMEO2","CHCL2CO3","CLCOCH2O2","CLCOCLO2","CCLOHCO3","HNMVKO2",
"NC3CO3","C42O2","HC3CO3","C41O2","MVKOHAO2","MVKOHBO2","HC3CCO3","INCO2","NC4CO3","C510O2","C536O2","C537O2","INAO2","C58O2","HC4CO3",
"CHOCOMOXO2","ACETMECO3","HOACETETO2","MECOACETO2","ACPRONEO2","ACCOETO2","ACETC2CO3","IPRACBCO3","ACBUONEAO2","ACBUONEBO2",
"ACCOC3H6O2","SBUACONEO2","TBUACCO3","MTBEACHOO2","MTBEBCO3","IPROC21O2","IPROMCCO3","EIPEO2","ETBEACO3","ETBECCO3","BOXCOEOLO2",
"BRETO3","HO1CO3C4O2","BIACETO2","HO2CO4C5O2","CO23C54O2","HOCO3C54O2","C53O2","C41CO3","CO2HOC61O2","CO24C6O2","CO25C6O2","HO2C4O2",
"C61O2","CO23C65O2","C6CO3OH5O2","C6CO34O2","C6HO1CO3O2","C3COCCO3","PEN2ONE1O2","MIBK3COO2","C612O2","CO2M33CO3","C6COCHOO2",
"CY6DIONO2","NC101O2","C96CO3","C720O2","NC91CO3","C8BCO2","C918CO3","C923CO3","C141CO3","NBCALO2","BCALO2","BCSOZO2","C151O2",
"C152O2","MMFO2","MMCFO2","DMSO2O2","CHOC4CO3","C6DIALO2","CHOC4O2","CYC6DIONO2","CONM2CO3","NBZFUO2","BZFUO2","CATEC1O2","MCATEC1O2",
"MC3CODBCO3","C4M2ALOHO2","C5DICARBO2","NTLFUO2","TLFUO2","MC4CODBCO3","MC5CO2OHO2","NOXYFUO2","C6OTKETO2","OXYFUO2","OXCATEC1O2",
"C5MCO2OHO2","NMXYFUO2","C23O3MO2","MXYFUO2","NPXYFUO2","MCOCOMOXO2","PXYFUO2","MXCATEC1O2","DMKOHO2","PXCATEC1O2","ECATEC1O2",
"C6DICARBO2","NEBFUO2","BUTALAO2","EBFUO2","C7CO3OHO2","PCATEC1O2","C7DCO2","NPBFUO2","C4CHOAO2","PBFUO2","C8CO3OHO2","PHCOMEO2",
"IPCATEC1O2","IC7DCO2","NIPBFUO2","IC4CHOAO2","IPBFUO2","C7MCO3OHO2","T123CAT1O2","C7ADCCO3","C7ADCO2","NTMB1FUO2","TMB1FUO2",
"NTMB2FUO2","MC6OTKETO2","TMB2FUO2","C7BDCO2","T124CAT1O2","OTCATEC1O2","MTCATEC1O2","C7EDCO2","PTCATEC1O2","C7DDCCO3","C7DDCO2",
"NMEBFUO2","C23O3EO2","MEBFUO2","EMPHO2","CH3COCCLO2","CLCOCCL2O2","C527O2","C526O2","HC4ACO3","C58AO2","INB1O2","INB2O2","HPC52O2",
"HC4CCO3","C57AO2","C57O2","INDO2","C59O2","C524O2","ETHFORMO2","IPRMEETO2","CHOOMCO3","PRONFORMO2","PRCOOMCO3","PRCOOMO2",
"BOXCOCHOO2","BOXFORMO2","PRONOCOPO2","BOXCOCOMO2","PINALO2","C108O2","C89CO3","C920CO3","C920O2","C97O2","C85CO3","C85O2","C719O2",
"C918O2","C9DCO2","C915O2","C917O2","NLIMALO2","LIMALO2","C729CO3","C822CO3","C924O2","C816CO3","NORLIMO2","C816O2","NLMKAO2",
"LMKAO2","LMKBO2","C146O2","C131CO3","BCLKAO2","BCLKBO2","BCLKCO2","C131O2","C147O2","C126CO3","C136CO3","C148O2","C1311O2",
"NC1313O2","C1313O2","C126O2","C144O2","C142O2","NBCKO2","BCKAO2","BCKBO2","CH3SOO","H13C43CO3","C42CO3","HOC3H6CO3","C3DIOLO2",
"HO2C4CO3","HOIBUTCO3","C63O2","HO3C5CO3","C64O2","HO2C54O2","HO2C54CO3","C66O2","CO3C54CO3","H2M2C4CO3","C67O2","C610O2","H2M3C4CO3",
"C68O2","C69O2","C611O2","HM33C3CO3","HM22C3O2","HM22C3CO3","HM2C43CO3","C71O2","C76O2","C77O2","C78O2","C711O2","H3M3C5O2",
"H3M3C5CO3","C82O2","C81O2","C93O2","C92O2","HO6C7O2","C103O2","C102O2","HO7C8O2","C113O2","C112O2","HO8C9O2","C123O2","C122O2",
"CO1H63O2","C3NO3COO2","NPHEN1O2","NNCATECO2","NCATECO2","NBZQO2","PBZQO2","NPTLQO2","PTLQO2","NCRES1O2","MNNCATECO2","MNCATECO2",
"NOXYOL1O2","NOXYQO2","OXYQO2","OXNNCATCO2","OXNCATECO2","C534O2","NMXYOL1O2","NMXYQO2","MXYQO2","MXNNCATCO2","MXNCATECO2",
"NPXYOL1O2","NPXYQO2","PXYQO2","PXNNCATCO2","PXNCATECO2","NEBNZ1O2","NPEBQO2","PEBQO2","ENNCATECO2","ENCATECO2","CO3H4CO3",
"PHCOCOCO2","NPBNZ1O2","NPPRBQO2","PPRBQO2","PNNCATECO2","PNCATECO2","C5O45OHCO3","NIPBNZ1O2","NIPRBQO2","IPRBQO2","IPNNCATCO2",
"IPNCATECO2","C4MOHOCO3","NT123L1O2","T123NNCTO2","T123NCATO2","NT124L1O2","NTM124QO2","TM124QO2","T124NNCTO2","T124NCATO2",
"C5CO234O2","NOETOL1O2","NOETLQO2","OETLQO2","OTNNCATCO2","OTNCATECO2","NMETOL1O2","NMETLQO2","METLQO2","MTNNCATCO2","MTNCATECO2",
"NPETOL1O2","NPETLQO2","PETLQO2","PTNNCATCO2","PTNCATECO2","CO234C65O2","H13CO2CO3","CO2N3CO3","C535O2","C58NO3CO3","ACCOCOMEO2",
"ACEETOHO2","ACCOMCOMO2","ACCOCOETO2","MTBEAALCO3","C62O2","HO13C4O2","HM22CO3","C6COCHOCO3","C5COCHOO2","CHOC2H4O2","HCOCH2CO3",
"CY6TRIONO2","C6CYTONO2","NC102O2","C512CO3","C89O2","C926O2","C817CO3","C817O2","NC826O2","C826O2","C729O2","LMLKAO2","LMLKBO2",
"C116CO3","C116O2","C129O2","C1210O2","CH3SOO2","C1H4C5CO3","CHOC4OHO2","HOC4CHOO2","C6145COO2","COHM2CO3","CO2C4CO3","HOBUT2CO3",
"CO3C5CO3","CO2C54CO3","C65O2","CO2M3C4CO3","C72O2","CO25C73O2","CO25C74O2","C712O2","C713O2","C714O2","C84O2","C94O2","C104O2",
"C114O2","C6H13CO3","C124O2","MALANHYO2","NDNPHENO2","DNPHENO2","NDNCRESO2","DNCRESO2","C6O4KETO2","NDNOXYOLO2","DNOXYOLO2",
"MMALANHYO2","CH3COCO3","NDNMXYOLO2","DNMXYOLO2","TL4OHNO2O2","NDNPXYOLO2","DNPXYOLO2","NDNEBNZLO2","DNEBNZLO2","NDNPBNZLO2",
"DNPBNZLO2","C61CO3","NDNIPBZLO2","DNIPBNZLO2","C62CO3","NDNT123LO2","DNT123LO2","TM124NO2O2","NDNT124LO2","DNT124LO2","MXYOHNO2O2",
"NDNOETOLO2","DNOETOLO2","NDNMETOLO2","DNMETOLO2","NDNPETOLO2","DNPETOLO2","CO356OCO2","C531O2","INCNCO3","IEACO3","IECCO3","HPC52CO3",
"INDHCO3","C57NO3CO3","INAHPCO3","INANCO3","INAHCO3","NC524O2","C525O2","HMACO3","HMACRO2","ACCOMECO3","IPRFORMO2","PRCOFORMO2",
"PRONOCOMO2","CO23C4CO3","C5CO34CO3","C106O2","C717O2","C811CO3","C921O2","C98O2","C86O2","C919O2","C914O2","C916O2","C88CO3",
"C88O2","C512O2","C619O2","C626CO3","C626O2","C735O2","C822O2","C823CO3","C925O2","C622CO3","C1011CO3","C1210CO3","C132O2","C137CO3",
"C1013CO3","C1312O2","C127O2","C143O2","CH3SO2O2","HO24C5O2","C55O2","C67CO3","H3M2C4CO3","C79O2","H3M3C4CO3","H13M3C5O2","HO4C5CO3",
"HO5C6CO3","HO6C7CO3","HO7C8CO3","HO8C9CO3","C5CO2OHCO3","C6CO2OHCO3","C5M2OHOCO3","C4COMOHCO3","C23O3MCO3","C23O3CCO3","C7CO2OHCO3",
"C6MOHCOCO3","C7OHCO2CO3","ECO3CO3","C8OHCO2CO3","C8CO2OHCO3","NDMMALYO2","DMMALYO2","C7MOHCOCO3","C5MEJCO3","C6EO2OHCO3","C7MJPCO3",
"C23O3ECO3","EMPOHNO2O2","C47CO3","INB1HPCO3","INB1NACO3","INB1NBCO3","MMALNACO3","MMALNBCO3","INDHPCO3","INANCOCO3","HIEB1O2",
"HIEB2O2","HO13C3CO3","C5CO23O2","CHOC2CO3","CHOC3COCO3","C5124COCO3","CO235C6CO3","NC71O2","C811O2","CHOC3COO2","H3C25C6CO3",
"H3C25C6O2","C810O2","C818O2","C727CO3","NC728O2","C728O2","C622O2","C823O2","C819O2","C731CO3","C1011O2","C137O2","C1013O2",
"C1010O2","C117O2","C830CO3","C145O2","C927O2","C1214O2","CHOC4DOLO2","C6TRONOHO2","C23C54CO3","C73O2","C74O2","C715O2","C83O2",
"C95O2","C105O2","C115O2","C125O2","C4CO2DBCO3","C5CO2DBCO3","C4DBM2CO3","C5DBCO2CO3","C7CO2DBCO3","C8CO2DBCO3","C8DBCO2CO3",
"C4DBMECO3","C5DBECO3","C5EDBCO3","C31CO3","C533O2","MECOFORMO2","C5124COO2","CO235C6O2","C716O2","C922O2","C614O2","C511O2",
"C620O2","C87CO3","C616O2","C718CO3","C513O2","CO25C6CO3","C627O2","C727O2","C511CO3","C517CO3","C517O2","C628O2","C824O2","C1211CO3",
"C133O2","C830O2","C128O2","HO24C4CO3","C613O2","CO2OH3MCO3","C812O2","C721CO3","C721O2","H3C2C4CO3","C87O2","C718O2","C514O2",
"C820O2","C518CO3","NC623O2","C623O2","C825O2","C731O2","C732CO3","C1012O2","C1211O2","C139O2","C1014O2","C736O2","C118O2",
"C928CO3","C630O2","C1215O2","EMALANHYO2","PMALANHYO2","IPMALNHYO2","C312COCO3","CHOCOCH2O2","NC72O2","C621O2","C515CO3","C515O2",
"C821O2","HMVKBCO3","C520O2","C624CO3","C732O2","C829O2","C134O2","C827CO3","C522CO3","C831O2","C928O2","C46CO3","C930O2","C813O2",
"C722O2","C615CO3","C617CO3","C618CO3","C617O2","C618O2","NC730O2","C730O2","C624O2","C733O2","C1212O2","C827O2","C1310O2","NC61CO3",
"C615O2","C519CO3","C629O2","C734O2","C521O2","C135O2","COO2C4CO3","COO2C4O2","C929O2","C516O2","C44O2","H1C23C4CO3","H1C23C4O2",
"CO1M22CO3","C519O2","C625O2","C1213O2","COO2C3CO3","C828CO3","C828O2"]

function stoich_str2int(stoich::String)::Int
    if stoich==""
        return 1
    else
        return parse(Int,stoich)
    end
end

function parse_reactants(file::String)#,RO2_names::Array{String,1}
    reactants_dict=Dict{String,Set{Tuple{Int,Int}}}() # Dict{Reactants,Set{Tuple{LineNum,Stoich}}}
    reactants2ind=Dict{String,Int}()
    ind2reactants=Dict{Int,String}()
    reactant_step=0
    #file="MCM_test.eqn.txt"
    num_eqns=0
    open(file,"r") do f
        #global num_eqns
        for line in eachline(f)
            ind = parse(Int,match(r"\{([0-9]+)\.\}",line,1).captures[1])
            num_eqns=max(ind,num_eqns)
            line_parts=split(strip(line)[1:end-1],":")
            eqn=strip(split(line_parts[1],"}")[2])
            lhs_rhs=[strip(i) for i in split(eqn,"=")]
            lhs=[strip(i) for i in split(lhs_rhs[1],"+")]
            rhs=[strip(i) for i in split(lhs_rhs[2],"+")]
            reac_regexp=r"([0-9]*)([A-Z][0-9A-Z]*)"
            reacs=[Tuple{String,String}(match(reac_regexp,i).captures) for i in lhs]
            prods=[Tuple{String,String}(match(reac_regexp,i).captures) for i in rhs]
            for (stoich,reac) in reacs
                num_stoich=stoich_str2int(stoich)
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
            for (stoich,reac) in prods
                num_stoich=stoich_str2int(stoich)*(-1)# because it is product
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
            #println(reacs)
            #println(prods)
        end
    end

    num_reactants=length(reactants_dict)
    assert(num_reactants==reactant_step)
    #reactants_list=[ for i in 1:num_reactants]
    reactants_inds=1:num_reactants
    #reactants2ind=Dict(reactants_list[i]=>i for i in reactants_inds)

    RO2_inds=[reactants2ind[i] for i in RO2_names if haskey(reactants_dict,i)]

    stoich_mtx=spzeros(num_reactants,num_eqns)
    for i in reactants_inds
        for (eqn_ind,stoich) in reactants_dict[ind2reactants[i]]
            stoich_mtx[i,eqn_ind]=stoich
        end
    end
    dropzeros!(stoich_mtx)
    return (stoich_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind)
end

function gen_evaluate_rates(file)
    rate_expr=quote end
    #file="MCM_test.eqn.txt"
    open(file,"r") do f
        for line in eachline(f)
            line_parts=split(strip(line)[1:end-1],":")
            spd=strip(line_parts[2])
            #replacing 123.456D-12 to 123.456E-12
            num_regexp=r"([0-9]*\.?[0-9]+)D([\+\-]?[0-9]+)"
            spd_exprs=replace(spd,num_regexp,s"\1E\2")
            spd_exprs=replace(spd_exprs,"TEMP","temp")
            spd_exprs=replace(spd_exprs,"EXP","exp")
            spd_exprs=replace(spd_exprs,"**","^")
            spd_exprs=replace(spd_exprs,r"J\(([0-9]+)\)",s"J[\1]")
            #println(spd_expr)
            spd_expr=parse(spd_exprs)
            ind = parse(Int,match(r"\{([0-9]+)\.\}",line,1).captures[1])
            push!(rate_expr.args,:(rate_values[$ind]=$spd_expr))
        end
    end
    #op=parse(expr)
    return quote


function evaluate_rates!(ttime::Float64,RO2::Float64,H2O::Float64,temp::Float64,rate_values::Array{Float64,1},J::Array{Float64,1})
    # ttime: seconds
    # RO2: sum of RO2s
    # Now cycle through all of the mcm_constants_dict values.
    # Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    # need to modify this potentially.
    # Creating reference to constant values used in rate expressions 

    # constants used in calculation of reaction rates 
    M  = 2.55E+19  #Check this against pressure assumed in model 
    N2 = 0.79*M 
    O2 = 0.2095*M 

    # kro2no : ro2      + no      = ro      + no2 
    #        : ro2      + no      = rono2 
    # iupac 1992 
    KRONO2    = 2.70E-12*exp(360.0/temp) 
    
    KRO2NO = 2.7E-12*exp(360.0/temp)  

    # kro2ho2: ro2      + ho2     = rooh    + o2 
    # mcm protocol [1997] 
    KRO2HO2   = 2.91E-13*exp(1300.0/temp) 

    # kapho2 : rcoo2    + ho2     = products 
    # mcm protocol [1997] 
    KAPHO2    = 5.2E-13*exp(980.0/temp) 

    # kapno  : rcoo2    + no      = products 
    # mej [1998] 
    KAPNO = 7.5E-12*exp(290.0/temp) 

    # kro2no3: ro2      + no3     = products 
    # mcm protocol [1997] 
    KRO2NO3   = 2.3E-12 

    # kno3al : no3      + rcho    = rcoo2   + hno3 
    # mcm protocol [1997] 
    KNO3AL    = 1.4E-12*exp(-1860.0/temp) 

    # kdec   : ro                 = products 
    # mcm protocol [1997] 
    KDEC      = 1.00E+06 

    KROSEC = 2.50E-14*exp(-300.0/temp) 

    KALKOXY=3.70E-14*exp(-460.0/temp)*O2 

    KALKPXY=1.80E-14*exp(-260.0/temp)*O2 

    KROPRIM = 2.50E-14*exp(-300.0/temp) 

    KCH3O2 = 1.03E-13*exp(365.0/temp) 

    K298CH3O2 = 3.5E-13 

    K14ISOM1 = 3.00E7*exp(-5300.0/temp)

    # ------------------------------------------------------------------- 
    # complex reactions 
    # ------------------------------------------------------------------- 

    # kfpan kbpan 
    # formation and decomposition of pan 
    # iupac 2001 (mcmv3.2) 
    kc0     = 3.28E-28*M*(temp/300.0)^(-6.87) 
    kci     = 1.125E-11*(temp/300.0)^(-1.105) 
    krc     = kc0/kci 
    fcc     = 0.30 
    nc      = 0.75-(1.27*log10(fcc)) 
    fc      = 10^(log10(fcc)/(1.0+((log10(krc))/nc)^2.0)) 
    KFPAN   = (kc0*kci)*fc/(kc0+kci) 

    kd0     = 1.10E-05*M*exp(-10100.0/temp) 
    kdi     = 1.90E+17*exp(-14100.0/temp) 
    krd     = kd0/kdi 
    fcd     = 0.30 
    ncd     = 0.75-(1.27*log10(fcd)) 
    fd      = 10.0^(log10(fcd)/(1.0+((log10(krd))/ncd)^2.0)) 
    KBPAN   = (kd0*kdi)*fd/(kd0+kdi) 

    KPPN0     = 1.7E-03*M*exp(-11280.0/temp) 
    KPPNI     = 8.3E+16*exp(-13940.0/temp) 
    KRPPN     = KPPN0/KPPNI 
    FCPPN     = 0.36 
    NCPPN     = 0.75-(1.27*log10(fcd)) 
    FPPN      = 10.0^(log10(fcd)/(1.0+((log10(krd))/ncd)^2.0)) 
    KBPPN   = (KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI) 

    # kmt01  : o        + no      = no2 
    # iupac 2001 (mcmv3.2) 
    k10     = 1.00E-31*M*(temp/300.0)^(-1.6) 

    k1i     = 5.0E-11*(temp/300.0)^(0.3) 
    kr1     = k10/k1i 
    fc1     = 0.85 
    nc1     = 0.75-(1.27*log10(fc1)) 
    f1      = 10.0^(log10(fc1)/(1.0+((log10(kr1)/nc1))^2.0)) 
    KMT01   = (k10*k1i)*f1/(k10+k1i) 

    # kmt02  : o        + no2     = no3 
    # iupac 2001 (mcmv3.2) 
    k20     = 1.30E-31*M*(temp/300.0)^(-1.5) 
    k2i     = 2.30E-11*(temp/300.0)^(0.24) 
    kr2     = k20/k2i 
    fc2     = 0.6 
    nc2     = 0.75-(1.27*log10(fc2)) 
    f2      = 10.0^(log10(fc2)/(1.0+((log10(kr2)/nc2))^2.0)) 
    KMT02   = (k20*k2i)*f2/(k20+k2i) 

    # kmt03  : no2      + no3     = n2o5 
    # iupac 2006, mcmv3.2 
    k30     = 3.60E-30*M*(temp/300.0)^(-4.1) 
    k3i     = 1.90E-12*(temp/300.0)^(0.2) 
    kr3     = k30/k3i 
    fc3     = 0.35 
    nc3     = 0.75-(1.27*log10(fc3)) 
    f3      = 10.0^(log10(fc3)/(1.0+((log10(kr3)/nc3))^2.0)) 
    KMT03   = (k30*k3i)*f3/(k30+k3i) 

    # kmt04  : n2o5               = no2     + no3 
    # iupac 2006, mcmv3.2 
    k40     = 1.30E-03*M*(temp/300.0)^(-3.5)*exp(-11000.0/temp) 
    k4i     = 9.70E+14*(temp/300.0)^(0.1)*exp(-11080.0/temp) 
    kr4     = k40/k4i 
    fc4     = 0.35 
    nc4     = 0.75-(1.27*log10(fc4)) 
    f4      = 10.0^(log10(fc4)/(1+((log10(kr4)/nc4))^2.0)) 
    KMT04   = (k40*k4i)*f4/(k40+k4i) 

    # kmt05  : oh       + co(+o2) = ho2     + co2 
    # iupac 2006 
    KMT05  = 1.44E-13*(1.0 + (M/4.2E19)) 

    # kmt06  : ho2      + ho2     = h2o2    + o2 
    # water enhancement factor 
    # iupac 1992 
    KMT06  = 1.0 + (1.40E-21*exp(2200.0/temp)*H2O) 

    # kmt06  = 1.0 + (2.00E-25*exp(4670.0/temp)*h2o) 
    # S+R 2005 values 

    # kmt07  : oh       + no      = hono 

    # iupac 2006, mcmv3.2 
    k70     = 7.40E-31*M*(temp/300.0)^(-2.4) 
    k7i     = 3.30E-11*(temp/300.0)^(-0.3) 
    kr7     = k70/k7i 
    fc7     = 0.81 
    nc7     = 0.75-(1.27*log10(fc7)) 
    f7      = 10.0^(log10(fc7)/(1+((log10(kr7)/nc7))^2.0)) 
    KMT07   = (k70*k7i)*f7/(k70+k7i) 

    # kmt08  : oh       + no2     = hno3 

    # iupac 2006, mcmv3.2 
    k80     = 3.2E-30*M*(temp/300.0)^(-4.5) 
    k8i     = 3.0E-11 
    kr8     = k80/k8i 
    fc8     = 0.41 
    nc8     = 0.75-(1.27*log10(fc8)) 
    f8      = 10.0^(log10(fc8)/(1.0+((log10(kr8)/nc8))^2.0)) 
    KMT08   = (k80*k8i)*f8/(k80+k8i) 

    # kmt09  : ho2      + no2     = ho2no2 
    # iupac 1997, mcmv3.2 
    
    k90     = 1.4E-31*M*(temp/300.0)^(-3.1) 
    k9i     = 4.0E-12 
    kr9     = k90/k9i 
    fc9     = 0.4 
    nc9     = 0.75-(1.27*log10(fc9)) 
    f9      = 10.0^(log10(fc9)/(1.0+((log10(kr9)/nc9))^2.0)) 
    KMT09   = (k90*k9i)*f9/(k90+k9i) 

    # kmt10  : ho2no2             = ho2     + no2 
    # iupac 1997, mcmv3.2 

    k100     = 4.10E-05*M*exp(-10650.0/temp) 
    k10i     = 6.0E+15*exp(-11170.0/temp) 
    kr10     = k100/k10i 
    fc10     = 0.4 
    nc10     = 0.75-(1.27*log10(fc10)) 
    f10      = 10.0^(log10(fc10)/(1.0+((log10(kr10)/nc10))^2.0)) 
    KMT10    = (k100*k10i)*f10/(k100+k10i) 

    # kmt11  : oh       + hno3    = h2o     + no3 
    # iupac 2006, mcmv3.2 

    k1     = 2.40E-14*exp(460.0/temp) 
    k3     = 6.50E-34*exp(1335.0/temp) 
    k4     = 2.70E-17*exp(2199.0/temp) 
    k2     = (k3*M)/(1.0+(k3*M/k4)) 
    KMT11  = k1 + k2 

    # kmt12 iupac 2006, mcmv3.2 

    k120 = 2.5E-31*((temp/300.0)^(-2.6))*M 
    k12i = 2.0E-12 
    kr12 = k120/k12i 
    fc12 = 0.53 
    nc12 = 0.75-(1.27*log10(fc12)) 
    f12  = 10.0^(log10(fc12)/(1.0+((log10(kr12)/nc12))^2.0)) 
    KMT12    = (k120*k12i)*f12/(k120+k12i) 

    # kmt13  : ch3o2    + no2     = ch3o2no2 
    # iupac 2006 

    k130     = 2.50E-30*((temp/300.0)^(-5.5))*M 
    k13i     = 1.80E-11 
    kr13     = k130/k13i 
    fc13     = 0.36 
    nc13     = 0.75-(1.27*log10(fc13)) 
    f13      = 10.0^(log10(fc13)/(1.0+((log10(kr13)/nc13))^2.0)) 
    KMT13    = (k130*k13i)*f13/(k130+k13i) 

    # kmt14  : ch3o2no2           = ch3o2   + no2 
    # iupac 2006, mcmv3.2 

    k140     = 9.00E-05*exp(-9690.0/temp)*M 
    k14i     = 1.10E+16*exp(-10560.0/temp) 
    kr14     = k140/k14i 
    fc14     = 0.36 
    nc14     = 0.75-(1.27*log10(fc14)) 
    f14      = 10.0^(log10(fc14)/(1.0+((log10(kr14)/nc14))^2.0)) 
    KMT14    = (k140*k14i)*f14/(k140+k14i) 

    # kmt15 iupac 2006, mcmv3.2 

    k150 = 8.60E-29*((temp/300.0)^(-3.1))*M 
    k15i = 9.00E-12*((temp/300.0)^(-0.85)) 
    kr15 = k150/k15i 
    fc15 = 0.48 
    nc15 = 0.75-(1.27*log10(fc15)) 
    f15  = 10.0^(log10(fc15)/(1.0+((log10(kr15)/nc15))^2.0)) 
    KMT15 = (k150*k15i)*f15/(k150+k15i) 

    # kmt16  :  oh  +  c3h6 
    # iupac 2006 

    k160     = 8.00E-27*((temp/300.0)^(-3.5))*M 
    k16i     = 3.00E-11*((temp/300.0)^(-1.0)) 
    kr16     = k160/k16i 
    fc16     = 0.5 
    nc16     = 0.75-(1.27*log10(fc16)) 
    f16      = 10.0^(log10(fc16)/(1.0+((log10(kr16)/nc16))^2.0)) 
    KMT16    = (k160*k16i)*f16/(k160+k16i) 

    # kmt17 iupac 2006 

    k170 = 5.00E-30*((temp/300.0)^(-1.5))*M 
    k17i = 1.00E-12 
    kr17 = k170/k17i 
    fc17 = (0.17*exp(-51./temp))+exp(-1.0*temp/204.0) 
    nc17 = 0.75-(1.27*log10(fc17)) 
    f17  = 10.0^(log10(fc17)/(1.0+((log10(kr17)/nc17))^2.0)) 
    KMT17 = (k170*k17i)*f17/(k170+k17i) 
    
    KMT18 = 9.5E-39*O2*exp(5270.0/temp)/(1+7.5E-29*O2*exp(5610.0/temp)) 

    # ************************************************************************ 
    # define photolysis reaction rates using derwent method from mcm2box.fac 
    # ************************************************************************ 

    # solar declination angle  
    dec = 23.79 
    # latitude 
    lat = 50.0 
    pi = 4.0*atan(1.0) 
    # local hour angle - representing time of day 
    lha = (1.0+ttime/4.32E+4)*pi 
    radian = 180.0/pi 
    lat = lat/radian 
    dec = dec/radian 
    theta = acos(cos(lha)*cos(dec)*cos(lat)+sin(dec)*sin(lat)) 
    sinld = sin(lat)*sin(dec) 
    cosld = cos(lat)*cos(dec) 
    cosx = (cos(lha)*cosld)+sinld 
    cosx = cos(theta) 
    secx = 1.0E+0/(cosx+1.0E-30) 

    # Data taken from photolysis.txt. Calculations done in the form of: 
    # j(k) = l(k)*cosx**( mm(k))*exp(-nn(k)*secx) 
    #J=[None]*62 
    #J          L           M          N 
    J[1]=6.073E-05*cosx^(1.743)*exp(-1.0*0.474*secx) 
    J[2]=4.775E-04*cosx^(0.298)*exp(-1.0*0.080*secx) 
    J[3]=1.041E-05*cosx^(0.723)*exp(-1.0*0.279*secx) 
    J[4]=1.165E-02*cosx^(0.244)*exp(-1.0*0.267*secx) 
    J[5]=2.485E-02*cosx^(0.168)*exp(-1.0*0.108*secx) 
    J[6]=1.747E-01*cosx^(0.155)*exp(-1.0*0.125*secx) 
    J[7]=2.644E-03*cosx^(0.261)*exp(-1.0*0.288*secx) 
    J[8]=9.312E-07*cosx^(1.230)*exp(-1.0*0.307*secx) 
    J[11]=4.642E-05*cosx^(0.762)*exp(-1.0*0.353*secx) 
    J[12]=6.853E-05*cosx^(0.477)*exp(-1.0*0.323*secx) 
    J[13]=7.344E-06*cosx^(1.202)*exp(-1.0*0.417*secx) 
    J[14]=2.879E-05*cosx^(1.067)*exp(-1.0*0.358*secx) 
    J[15]=2.792E-05*cosx^(0.805)*exp(-1.0*0.338*secx) 
    J[16]=1.675E-05*cosx^(0.805)*exp(-1.0*0.338*secx) 
    J[17]=7.914E-05*cosx^(0.764)*exp(-1.0*0.364*secx) 
    J[18]=1.482E-06*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[19]=1.482E-05*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[20]=7.600E-04*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[21]=7.992E-07*cosx^(1.578)*exp(-1.0*0.271*secx) 
    J[22]=5.804E-06*cosx^(1.092)*exp(-1.0*0.377*secx) 
    J[23]=2.4246E-06*cosx^(0.395)*exp(-1.0*0.296*secx) 
    J[24]=2.424E-06*cosx^(0.395)*exp(-1.0*0.296*secx) 
    J[31]=6.845E-05*cosx^(0.130)*exp(-1.0*0.201*secx) 
    J[32]=1.032E-05*cosx^(0.130)*exp(-1.0*0.201*secx) 
    J[33]=3.802E-05*cosx^(0.644)*exp(-1.0*0.312*secx) 
    J[34]=1.537E-04*cosx^(0.170)*exp(-1.0*0.208*secx) 
    J[35]=3.326E-04*cosx^(0.148)*exp(-1.0*0.215*secx) 
    J[41]=7.649E-06*cosx^(0.682)*exp(-1.0*0.279*secx) 
    J[51]=1.588E-06*cosx^(1.154)*exp(-1.0*0.318*secx) 
    J[52]=1.907E-06*cosx^(1.244)*exp(-1.0*0.335*secx) 
    J[53]=2.485E-06*cosx^(1.196)*exp(-1.0*0.328*secx) 
    J[54]=4.095E-06*cosx^(1.111)*exp(-1.0*0.316*secx) 
    J[55]=1.135E-05*cosx^(0.974)*exp(-1.0*0.309*secx) 
    J[56]=4.365E-05*cosx^(1.089)*exp(-1.0*0.323*secx) 
    J[57]=3.363E-06*cosx^(1.296)*exp(-1.0*0.322*secx) 
    J[61]=7.537E-04*cosx^(0.499)*exp(-1.0*0.266*secx) 
    
    $rate_expr
    return rate_values
end
    end
end