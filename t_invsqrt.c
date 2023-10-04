/*WBL 4 Oct 2023 use sensible name for table*/
/*GI produced data for log2*/
/*Created by make_t_log2.awk $Revision: 1.1 $ 10 Mar 2019 from gi_t_invsqrt.out8*/
/*gi_t_invsqrt.bat seed 140353 $Revision: 1.15 $ /cs/research/crest/projects1/ucacbbl/gi_glibc/invsqrt/cmaes_invsqrt eden.cs.ucl.ac.uk /tmp/gi_t_invsqrt Sun 10 Mar 14:41:11 GMT 2019*/

const float __t_invsqrt[512] = {
1.42094767093658447, 1.40222084522247314, 1.39880394935607910, 1.40915036201477051,
1.39340054988861084, 1.39251077175140381, 1.39951300621032715, 1.39333713054656982,
1.38411426544189453, 1.39703845977783203, 1.38343918323516846, 1.38024663925170898,
1.37570428848266602, 1.37972724437713623, 1.36949026584625244, 1.37159740924835205,
1.37128210067749023, 1.36747384071350098, 1.37148475646972656, 1.36498391628265381,
1.35768175125122070, 1.35909664630889893, 1.34840440750122070, 1.35890555381774902,
1.35664403438568115, 1.34027314186096191, 1.34396517276763916, 1.34593284130096436,
1.33906626701354980, 1.33458256721496582, 1.33544564247131348, 1.32569086551666260,
1.32334160804748535, 1.33238649368286133, 1.32166004180908203, 1.32839405536651611,
1.32998633384704590, 1.32756614685058594, 1.32579994201660156, 1.31084883213043213,
1.31197619438171387, 1.31923580169677734, 1.30899894237518311, 1.29968023300170898,
1.31199395656585693, 1.30041813850402832, 1.29894554615020752, 1.30703294277191162,
1.29286479949951172, 1.29657447338104248, 1.29162824153900146, 1.29629063606262207,
1.28400731086730957, 1.29389226436614990, 1.27768385410308838, 1.27511811256408691,
1.28364086151123047, 1.27862250804901123, 1.28450393676757812, 1.28144204616546631,
1.27955496311187744, 1.26135802268981934, 1.27128732204437256, 1.25693440437316895,
1.26417875289916992, 1.25510227680206299, 1.26411449909210205, 1.25667989253997803,
1.24823737144470215, 1.25379490852355957, 1.25207388401031494, 1.25815141201019287,
1.25215518474578857, 1.25160968303680420, 1.24582183361053467, 1.24552690982818604,
1.24121701717376709, 1.24419522285461426, 1.23545730113983154, 1.22824335098266602,
1.23431873321533203, 1.22812294960021973, 1.22342312335968018, 1.22860598564147949,
1.22296631336212158, 1.22259283065795898, 1.22399044036865234, 1.21345412731170654,
1.22530603408813477, 1.21978533267974854, 1.21650266647338867, 1.20693051815032959,
1.21215820312500000, 1.20480763912200928, 1.21437835693359375, 1.21502852439880371,
1.21275639533996582, 1.20875167846679688, 1.19550085067749023, 1.19325768947601318,
1.19358623027801514, 1.19217252731323242, 1.19927966594696045, 1.19657945632934570,
1.19870877265930176, 1.18929862976074219, 1.19037127494812012, 1.18478882312774658,
1.18618988990783691, 1.17970108985900879, 1.17812359333038330, 1.18108308315277100,
1.17556118965148926, 1.18169653415679932, 1.17902708053588867, 1.16756188869476318,
1.17160344123840332, 1.16401159763336182, 1.16996717453002930, 1.16737568378448486,
1.16586244106292725, 1.16061639785766602, 1.17129838466644287, 1.15898621082305908,
1.15835678577423096, 1.15380036830902100, 1.15251922607421875, 1.16231572628021240,
1.16207432746887207, 1.15387773513793945, 1.14722847938537598, 1.15632188320159912,
1.14588153362274170, 1.14112126827239990, 1.14557611942291260, 1.14710438251495361,
1.13966715335845947, 1.14326751232147217, 1.14547121524810791, 1.14320266246795654,
1.12892448902130127, 1.13682627677917480, 1.13340497016906738, 1.12676322460174561,
1.13348031044006348, 1.12960779666900635, 1.13062405586242676, 1.12487435340881348,
1.13169991970062256, 1.12728619575500488, 1.12313759326934814, 1.12456798553466797,
1.12627017498016357, 1.12638962268829346, 1.11031568050384521, 1.10806679725646973,
1.11468470096588135, 1.10764777660369873, 1.11485314369201660, 1.10641992092132568,
1.10741353034973145, 1.11300718784332275, 1.11024975776672363, 1.11016380786895752,
1.09996330738067627, 1.09618997573852539, 1.10393106937408447, 1.09322714805603027,
1.10084259510040283, 1.09330761432647705, 1.09101116657257080, 1.10151386260986328,
1.09639370441436768, 1.08435189723968506, 1.08544456958770752, 1.09558403491973877,
1.08562505245208740, 1.07891201972961426, 1.08365094661712646, 1.09153139591217041,
1.07807278633117676, 1.07956385612487793, 1.08743047714233398, 1.07279276847839355,
1.08416950702667236, 1.07887876033782959, 1.07865393161773682, 1.07862162590026855,
1.07279539108276367, 1.06742179393768311, 1.07248067855834961, 1.06566536426544189,
1.06091415882110596, 1.06911706924438477, 1.06826627254486084, 1.07101094722747803,
1.06309056282043457, 1.06936788558959961, 1.06896162033081055, 1.05710184574127197,
1.05308783054351807, 1.05914318561553955, 1.06406688690185547, 1.05271422863006592,
1.05421745777130127, 1.05986797809600830, 1.05822193622589111, 1.05597579479217529,
1.04637515544891357, 1.04856324195861816, 1.05117213726043701, 1.04569578170776367,
1.04465234279632568, 1.05022108554840088, 1.04169905185699463, 1.04678738117218018,
1.04326367378234863, 1.04021143913269043, 1.03986632823944092, 1.04149043560028076,
1.03848063945770264, 1.02899992465972900, 1.03727555274963379, 1.03998112678527832,
1.03329217433929443, 1.03348553180694580, 1.02367258071899414, 1.02957105636596680,
1.02109980583190918, 1.01912403106689453, 1.02505981922149658, 1.02203822135925293,
1.01753473281860352, 1.02447497844696045, 1.02526485919952393, 1.01673293113708496,
1.02088606357574463, 1.01602077484130859, 1.02325332164764404, 1.01745533943176270,
1.01894271373748779, 1.01634514331817627, 1.01571846008300781, 1.00790560245513916,
1.00570929050445557, 1.01238763332366943, 1.00567996501922607, 1.01130032539367676,
1.00579333305358887, 1.00601494312286377, 1.00220870971679688, 1.01066267490386963,
1.00951933860778809, 1.00728225708007812, 1.00756001472473145, 1.00598359107971191,
1.00063169002532959, 1.00179982185363770, 0.99668526649475098, 0.99656933546066284,
0.99551290273666382, 0.98973834514617920, 0.98672199249267578, 0.98254984617233276,
0.98366260528564453, 0.98805874586105347, 0.97541719675064087, 0.97635072469711304,
0.97233301401138306, 0.97119259834289551, 0.97781574726104736, 0.97601640224456787,
0.96796756982803345, 0.96872127056121826, 0.97248041629791260, 0.96860599517822266,
0.96783900260925293, 0.96603637933731079, 0.95922231674194336, 0.95132994651794434,
0.95840281248092651, 0.94836592674255371, 0.95312047004699707, 0.94983345270156860,
0.94351500272750854, 0.94776308536529541, 0.95031523704528809, 0.94023799896240234,
0.93911653757095337, 0.94460105895996094, 0.93328559398651123, 0.94328707456588745,
0.94120651483535767, 0.93226426839828491, 0.93811237812042236, 0.92922085523605347,
0.93492358922958374, 0.92700451612472534, 0.92996263504028320, 0.91902136802673340,
0.92386895418167114, 0.92148602008819580, 0.92388373613357544, 0.92203247547149658,
0.91833120584487915, 0.90909165143966675, 0.91123819351196289, 0.91005706787109375,
0.91655999422073364, 0.90778726339340210, 0.90716791152954102, 0.91255021095275879,
0.90260654687881470, 0.90948027372360229, 0.90045076608657837, 0.90265756845474243,
0.89966744184494019, 0.90029489994049072, 0.89376914501190186, 0.89124447107315063,
0.88997948169708252, 0.89822113513946533, 0.89273303747177124, 0.89535647630691528,
0.88404834270477295, 0.88819658756256104, 0.88071841001510620, 0.88097208738327026,
0.88814175128936768, 0.88365650177001953, 0.88207232952117920, 0.88248682022094727,
0.87356150150299072, 0.87194252014160156, 0.87966132164001465, 0.87055838108062744,
0.86878842115402222, 0.86712169647216797, 0.87510108947753906, 0.86523813009262085,
0.87023752927780151, 0.86613792181015015, 0.87004452943801880, 0.85755401849746704,
0.86152988672256470, 0.85971039533615112, 0.85978567600250244, 0.85907506942749023,
0.85467201471328735, 0.85119366645812988, 0.85628312826156616, 0.85835129022598267,
0.85018706321716309, 0.85453110933303833, 0.84906220436096191, 0.84809744358062744,
0.84740751981735229, 0.84401810169219971, 0.84832864999771118, 0.84280371665954590,
0.84130489826202393, 0.84728968143463135, 0.83779692649841309, 0.84437811374664307,
0.83292466402053833, 0.83809059858322144, 0.83591783046722412, 0.83332341909408569,
0.83344137668609619, 0.83263099193572998, 0.82901531457901001, 0.82775253057479858,
0.83341479301452637, 0.83270436525344849, 0.83126533031463623, 0.82802766561508179,
0.81924194097518921, 0.81968361139297485, 0.82815223932266235, 0.81611663103103638,
0.81650274991989136, 0.82407331466674805, 0.82280683517456055, 0.82023125886917114,
0.81429576873779297, 0.81225270032882690, 0.80909973382949829, 0.81662178039550781,
0.81397634744644165, 0.81378614902496338, 0.80686515569686890, 0.80892193317413330,
0.81085211038589478, 0.80918079614639282, 0.80053997039794922, 0.80496495962142944,
0.80328553915023804, 0.80173647403717041, 0.79789608716964722, 0.80440253019332886,
0.79597121477127075, 0.80147802829742432, 0.79953461885452271, 0.79758590459823608,
0.79586291313171387, 0.79204154014587402, 0.79616618156433105, 0.79531067609786987,
0.79317194223403931, 0.78990870714187622, 0.79450863599777222, 0.78451305627822876,
0.78510850667953491, 0.78652960062026978, 0.78938812017440796, 0.78152734041213989,
0.78367644548416138, 0.78085732460021973, 0.77920383214950562, 0.77682709693908691,
0.77769678831100464, 0.77944785356521606, 0.77858513593673706, 0.77621471881866455,
0.78100931644439697, 0.78032284975051880, 0.77008342742919922, 0.77697533369064331,
0.77519810199737549, 0.76735711097717285, 0.77114224433898926, 0.76527792215347290,
0.76457482576370239, 0.76925557851791382, 0.76626276969909668, 0.77118635177612305,
0.76835858821868896, 0.76681977510452271, 0.76856577396392822, 0.76342898607254028,
0.76535636186599731, 0.76205885410308838, 0.75903648138046265, 0.76087474822998047,
0.76008880138397217, 0.75296080112457275, 0.76015698909759521, 0.75910717248916626,
0.75640481710433960, 0.75159406661987305, 0.75128388404846191, 0.75606131553649902,
0.75625109672546387, 0.74728530645370483, 0.75284349918365479, 0.74510502815246582,
0.75270551443099976, 0.74834161996841431, 0.74327653646469116, 0.74924075603485107,
0.74570220708847046, 0.74616277217864990, 0.74382603168487549, 0.74359613656997681,
0.74201887845993042, 0.73794162273406982, 0.74376904964447021, 0.73849177360534668,
0.74105775356292725, 0.73945790529251099, 0.74174958467483521, 0.73368185758590698,
0.73692452907562256, 0.73912358283996582, 0.73120218515396118, 0.73587918281555176,
0.73165541887283325, 0.73529535531997681, 0.72879707813262939, 0.72727328538894653,
0.72885048389434814, 0.73033958673477173, 0.72910141944885254, 0.72854733467102051,
0.72633606195449829, 0.73072063922882080, 0.72690832614898682, 0.72135406732559204,
0.72619670629501343, 0.71820098161697388, 0.71922177076339722, 0.72540569305419922,
0.72270661592483521, 0.71631675958633423, 0.71804022789001465, 0.72142922878265381,
0.72238653898239136, 0.72260832786560059, 0.71911859512329102, 0.71754920482635498,
0.71946638822555542, 0.71690005064010620, 0.71453422307968140, 0.71031689643859863,
0.71135348081588745, 0.71115320920944214, 0.71325296163558960, 0.70653986930847168,
0.70714551210403442, 0.70564424991607666, 0.71214801073074341, 0.70562922954559326
};
/*gi_t_invsqrt.bat done Sun 10 Mar 14:41:17 GMT 2019*/
