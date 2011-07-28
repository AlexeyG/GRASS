#include "Configuration.h"
#include "ContigInfo.h"
#include "Reader.h"
#include "Writer.h"
#include <iostream>
#include <algorithm>

using namespace std;

Configuration config;
vector<ContigInfo> contigs;

bool readContigs(const string &fileName, vector<ContigInfo> &contigs)
{
	bool result = true;
	FastAReader reader;
	FastASequence seq;
	if (!reader.Open(fileName))
		result = false;
	while (reader.Read(seq))
		contigs.push_back(seq);
	reader.Close();
	sort(contigs.begin(), contigs.end());
	return result;
}

int getReadPosition(const string &str)
{
	int length = str.length();
	if (length == 0)
		return -1;
	int i = 0;
	while (i < length && str[i++] != '|');
	if (i >= length)
		return -1;
	int num = 0;
	for (; i < length; i++)
		num = num * 10 + str[i] - '0';
	return num;
}

bool inList(const string &name)
{
	string list = "16220.1|24871+16220.2|25377|148998.1|24924+148998.2|25430|235916.1|24639+235916.2|25142|298265.1|24927+298265.2|25426|227386.1|24884+227386.2|25383|194200.1|24765+194200.2|25259|190154.1|24657+190154.2|25148|375418.1|24726+375418.2|25217|229109.1|24785+229109.2|25276|363625.1|24948+363625.2|25438|104749.1|24938+104749.2|25427|138830.1|24674+138830.2|25162|313935.1|24852+313935.2|25340|298903.1|24867+298903.2|25354|21075.1|24916+21075.2|25401|214045.1|24672+214045.2|25157|235934.1|24723+235934.2|25206|33453.1|24985+33453.2|25468|74084.1|24663+74084.2|25146|81255.1|24670+81255.2|25153|130061.1|24799+130061.2|25282|367097.1|24651+367097.2|25133|149440.1|24841+149440.2|25323|278597.1|24774+278597.2|25256|325084.1|24710+325084.2|25192|130960.1|24723+130960.2|25204|243590.1|24670+243590.2|25151|122170.1|24806+122170.2|25287|138371.1|25017+138371.2|25498|221280.1|25010+221280.2|25491|371509.1|24895+371509.2|25376|82837.1|24761+82837.2|25241|183315.1|24971+183315.2|25451|76882.1|24675+76882.2|25155|337388.1|24631+337388.2|25111|14624.1|24658+14624.2|25137|366974.1|24773+366974.2|25252|322637.1|24902+322637.2|25381|79756.1|24732+79756.2|25211|102902.1|24879+102902.2|25357|52811.1|24680+52811.2|25158|137452.1|24815+137452.2|25293|37964.1|24930+37964.2|25408|276936.1|24980+276936.2|25457|215243.1|24934+215243.2|25411|268837.1|25008+268837.2|25485|345736.1|24647+345736.2|25124|261491.1|24993+261491.2|25470|22356.1|24888+22356.2|25365|140510.1|24903+140510.2|25380|123016.1|24840+123016.2|25317|237821.1|24660+237821.2|25137|29052.1|24725+29052.2|25202|188579.1|24649+188579.2|25125|332999.1|24660+332999.2|25136|175490.1|24788+175490.2|25264|307708.1|25019+307708.2|25494|164403.1|24942+164403.2|25417|103613.1|24736+103613.2|25211|121278.1|24998+121278.2|25473|271655.1|24753+271655.2|25228|36017.1|24767+36017.2|25241|145131.1|24676+145131.2|25150|40999.1|24808+40999.2|25282|246245.1|24729+246245.2|25203|89487.1|24838+89487.2|25312|141535.1|24870+141535.2|25344|279742.1|24801+279742.2|25275|117480.1|24851+117480.2|25325|94344.1|24881+94344.2|25354|237784.1|24882+237784.2|25355|168436.1|24777+168436.2|25250|43405.1|25017+43405.2|25490|278119.1|24658+278119.2|25131|246519.1|24989+246519.2|25462|156540.1|25013+156540.2|25485|217423.1|24954+217423.2|25426|106716.1|24981+106716.2|25453|343806.1|24775+343806.2|25247|260769.1|24683+260769.2|25155|246392.1|24790+246392.2|25262|97310.1|24866+97310.2|25338|39494.1|24758+39494.2|25229|320322.1|24988+320322.2|25459|208329.1|24755+208329.2|25226|228611.1|24921+228611.2|25392|212562.1|24762+212562.2|25233|239841.1|24655+239841.2|25125|305085.1|24737+305085.2|25207|286177.1|24711+286177.2|25181|145387.1|24713+145387.2|25183|150292.1|24658+150292.2|25128|81233.1|24743+81233.2|25213|180694.1|24660+180694.2|25130|322469.1|24879+322469.2|25349|374795.1|24813+374795.2|25282|208638.1|24691+208638.2|25160|135418.1|24706+135418.2|25175|87122.1|24653+87122.2|25122|295083.1|24687+295083.2|25156|388441.1|24757+388441.2|25225|9473.1|24731+9473.2|25199|245707.1|24835+245707.2|25303|355014.1|24674+355014.2|25142|36850.1|24783+36850.2|25251|106285.1|24749+106285.2|25217|59622.1|24954+59622.2|25422|292140.1|24871+292140.2|25339|57269.1|24670+57269.2|25137|99897.1|24692+99897.2|25159|55947.1|24972+55947.2|25439|338940.1|24764+338940.2|25231|42767.1|24798+42767.2|25265|121044.1|24664+121044.2|25131|39358.1|24730+39358.2|25197|27156.1|24865+27156.2|25332|368618.1|24937+368618.2|25404|187175.1|24952+187175.2|25419|316481.1|24698+316481.2|25165|36689.1|24854+36689.2|25320|249444.1|24941+249444.2|25407|244166.1|24884+244166.2|25350|361824.1|24712+361824.2|25178|349316.1|24717+349316.2|25183|115457.1|24851+115457.2|25317|270994.1|24860+270994.2|25326|107495.1|24903+107495.2|25369|231025.1|24736+231025.2|25202|168400.1|24928+168400.2|25394|335830.1|24793+335830.2|25259|201294.1|24987+201294.2|25453|211537.1|24677+211537.2|25143|163671.1|24763+163671.2|25228|376860.1|24898+376860.2|25363|233160.1|24685+233160.2|25150|341054.1|24995+341054.2|25460|27820.1|24968+27820.2|25433|193879.1|24882+193879.2|25347|31810.1|24904+31810.2|25369|165343.1|24836+165343.2|25301|254108.1|24924+254108.2|25389|46651.1|24970+46651.2|25434|112793.1|24702+112793.2|25166|270286.1|25002+270286.2|25466|116830.1|24997+116830.2|25461|194640.1|24837+194640.2|25301|313403.1|24675+313403.2|25139|29321.1|24856+29321.2|25320|139303.1|24761+139303.2|25225|158251.1|24781+158251.2|25245|186852.1|24688+186852.2|25152|102458.1|24662+102458.2|25126|280282.1|24865+280282.2|25329|105329.1|24859+105329.2|25323|107264.1|24953+107264.2|25417|2884.1|24766+2884.2|25229|301135.1|24897+301135.2|25360|140688.1|24670+140688.2|25133|199158.1|24898+199158.2|25361|138552.1|24891+138552.2|25354|337949.1|24865+337949.2|25328|240253.1|24664+240253.2|25127|219186.1|24656+219186.2|25119|133581.1|24910+133581.2|25373|332235.1|24749+332235.2|25212|45933.1|24742+45933.2|25205|160686.1|24685+160686.2|25148|329213.1|24765+329213.2|25227|336249.1|24981+336249.2|25443|91962.1|24908+91962.2|25370|332011.1|24860+332011.2|25322|311685.1|24921+311685.2|25383|305941.1|25010+305941.2|25472|135855.1|24768+135855.2|25230|318402.1|24878+318402.2|25339|146372.1|24715+146372.2|25176|12610.1|24827+12610.2|25288|77163.1|24738+77163.2|25199|194851.1|24748+194851.2|25209|104973.1|24895+104973.2|25356|287559.1|24989+287559.2|25450|209303.1|24912+209303.2|25373|69397.1|24897+69397.2|25358|215851.1|24726+215851.2|25187|89715.1|24699+89715.2|25160|38817.1|24748+38817.2|25208|190853.1|24769+190853.2|25229|369146.1|24994+369146.2|25454|235182.1|24990+235182.2|25450|211018.1|25013+211018.2|25473|235154.1|24918+235154.2|25378|278653.1|24824+278653.2|25284|290704.1|24949+290704.2|25409|221486.1|24784+221486.2|25244|303005.1|24902+303005.2|25361|372080.1|24795+372080.2|25254|98581.1|24702+98581.2|25161|110741.1|24803+110741.2|25262|21395.1|24892+21395.2|25351|366797.1|24805+366797.2|25264|277000.1|24717+277000.2|25176|238943.1|24839+238943.2|25298|70468.1|24768+70468.2|25227|133553.1|24725+133553.2|25184|163077.1|24708+163077.2|25167|351584.1|24813+351584.2|25271|78955.1|24934+78955.2|25392|103452.1|24759+103452.2|25217|62922.1|24705+62922.2|25163|86606.1|24917+86606.2|25375|264499.1|25023+264499.2|25481|74952.1|24921+74952.2|25379|173154.1|24833+173154.2|25290|214113.1|24666+214113.2|25123|217340.1|24864+217340.2|25321|186753.1|25014+186753.2|25471|196471.1|24924+196471.2|25381|97379.1|24686+97379.2|25143|139688.1|25019+139688.2|25476|111551.1|24840+111551.2|25297|87544.1|24752+87544.2|25209|365062.1|24654+365062.2|25111|40040.1|24975+40040.2|25432|120620.1|25018+120620.2|25475|98297.1|24655+98297.2|25111|174134.1|24988+174134.2|25444|278701.1|24987+278701.2|25443|321445.1|24816+321445.2|25272|198622.1|24678+198622.2|25134|155851.1|25011+155851.2|25467|46374.1|24816+46374.2|25271|370187.1|24660+370187.2|25115|344670.1|24719+344670.2|25174|294489.1|24785+294489.2|25240|61778.1|24850+61778.2|25305|98892.1|24889+98892.2|25344|105510.1|24956+105510.2|25411|277986.1|25005+277986.2|25460|174177.1|24768+174177.2|25223|366809.1|24916+366809.2|25371|31242.1|24857+31242.2|25312|241843.1|24769+241843.2|25224|156595.1|24762+156595.2|25217|26683.1|24998+26683.2|25452|24239.1|24704+24239.2|25158|22651.1|24669+22651.2|25123|328441.1|24774+328441.2|25228|321089.1|24672+321089.2|25126|314964.1|24750+314964.2|25204|9133.1|24711+9133.2|25165|90443.1|24828+90443.2|25282|101572.1|24989+101572.2|25443|215355.1|25014+215355.2|25468|256667.1|24785+256667.2|25239|8688.1|24956+8688.2|25410|120763.1|24698+120763.2|25152|130307.1|24844+130307.2|25298|154050.1|24680+154050.2|25134|377356.1|24830+377356.2|25283|357994.1|24968+357994.2|25421|356687.1|24785+356687.2|25238|313659.1|24782+313659.2|25235|216980.1|24865+216980.2|25318|223618.1|24895+223618.2|25348|160487.1|24896+160487.2|25349|146697.1|24724+146697.2|25177|140652.1|24980+140652.2|25433|106074.1|24926+106074.2|25379|103468.1|24697+103468.2|25150|102494.1|24685+102494.2|25138|74803.1|24986+74803.2|25439|72248.1|24685+72248.2|25138|62445.1|24990+62445.2|25443|171863.1|24771+171863.2|25224|56297.1|24658+56297.2|25111|182564.1|24946+182564.2|25399|208691.1|24717+208691.2|25169|267822.1|24757+267822.2|25209|376522.1|24818+376522.2|25270|104177.1|24902+104177.2|25354|21250.1|24879+21250.2|25331|25962.1|24966+25962.2|25418|379052.1|24661+379052.2|25113|386168.1|24857+386168.2|25309|195262.1|24975+195262.2|25427|56925.1|24841+56925.2|25293|337663.1|24898+337663.2|25350|222575.1|24708+222575.2|25160|345594.1|24704+345594.2|25156|137172.1|24802+137172.2|25253|240065.1|24806+240065.2|25257|253457.1|24734+253457.2|25185|261243.1|24778+261243.2|25229|386227.1|24998+386227.2|25449|116028.1|24823+116028.2|25274|368528.1|25020+368528.2|25471|95965.1|24939+95965.2|25390|91594.1|24960+91594.2|25411|303405.1|24733+303405.2|25184|40133.1|24890+40133.2|25341|350007.1|24792+350007.2|25243|348488.1|24836+348488.2|25287|55568.1|24662+55568.2|25113|186246.1|24983+186246.2|25434|183280.1|24710+183280.2|25161|99567.1|24761+99567.2|25211|361634.1|24918+361634.2|25368|385774.1|24702+385774.2|25152|110196.1|24786+110196.2|25236|144948.1|24676+144948.2|25126|28944.1|24884+28944.2|25334|153713.1|24924+153713.2|25374|115594.1|25020+115594.2|25470|209473.1|24813+209473.2|25263|247442.1|24715+247442.2|25165|71969.1|24815+71969.2|25265|357113.1|24969+357113.2|25418|342333.1|24748+342333.2|25197|66974.1|24890+66974.2|25339|373941.1|24970+373941.2|25419|261021.1|24683+261021.2|25132|267214.1|24773+267214.2|25222|272695.1|24981+272695.2|25430|215478.1|24977+215478.2|25426|836.1|24982+836.2|25431|127469.1|24998+127469.2|25446|245203.1|24950+245203.2|25398|165976.1|24785+165976.2|25233|102780.1|24664+102780.2|25112|195180.1|24967+195180.2|25414|128937.1|24898+128937.2|25345|87767.1|24761+87767.2|25208|343629.1|24950+343629.2|25397|48323.1|24976+48323.2|25423|189439.1|24752+189439.2|25199|379621.1|24725+379621.2|25172|212515.1|24830+212515.2|25277|14432.1|24792+14432.2|25239|97634.1|24777+97634.2|25223|294810.1|24965+294810.2|25411|355393.1|24838+355393.2|25284|93785.1|24979+93785.2|25425|185900.1|24742+185900.2|25188|167227.1|24798+167227.2|25244|334443.1|24720+334443.2|25166|173254.1|24783+173254.2|25228|194922.1|24898+194922.2|25343|246677.1|24944+246677.2|25389|157056.1|24694+157056.2|25139|220226.1|24810+220226.2|25255|165261.1|24694+165261.2|25139|305564.1|24683+305564.2|25128|171075.1|24858+171075.2|25303|33714.1|24767+33714.2|25212|140035.1|24729+140035.2|25173|261931.1|24957+261931.2|25401|160945.1|24769+160945.2|25213|352600.1|24952+352600.2|25396|249916.1|24975+249916.2|25419|200696.1|24698+200696.2|25142|162844.1|25000+162844.2|25444|367389.1|24690+367389.2|25134|10707.1|24699+10707.2|25142|153644.1|24683+153644.2|25126|167240.1|24809+167240.2|25252|148565.1|24946+148565.2|25389|351178.1|24806+351178.2|25249|325050.1|24772+325050.2|25215|260853.1|24860+260853.2|25303|167171.1|24870+167171.2|25313|112309.1|24837+112309.2|25280|262564.1|24908+262564.2|25351|293653.1|24934+293653.2|25376|311721.1|24963+311721.2|25405|270159.1|24810+270159.2|25252|6817.1|24884+6817.2|25326|274532.1|24762+274532.2|25204|275880.1|24683+275880.2|25125|278554.1|24681+278554.2|25123|160514.1|24712+160514.2|25154|79225.1|24672+79225.2|25114|134666.1|24937+134666.2|25379|159731.1|24976+159731.2|25418|65678.1|24813+65678.2|25255|154023.1|24674+154023.2|25115|243520.1|24825+243520.2|25266|283687.1|24686+283687.2|25127|248435.1|24714+248435.2|25155|144703.1|25016+144703.2|25457|85411.1|24828+85411.2|25269|351346.1|24793+351346.2|25234|188786.1|24997+188786.2|25438|78851.1|24724+78851.2|25165|375635.1|24826+375635.2|25266|148661.1|24944+148661.2|25384|168406.1|24836+168406.2|25276|19489.1|24812+19489.2|25252|214457.1|24765+214457.2|25205|302348.1|24831+302348.2|25271|298364.1|24901+298364.2|25341|96194.1|24860+96194.2|25299|340269.1|24682+340269.2|25121|200147.1|24784+200147.2|25223|94257.1|24949+94257.2|25388|353696.1|24961+353696.2|25400|64906.1|25001+64906.2|25440|67367.1|24825+67367.2|25264|135768.1|24893+135768.2|25332|130357.1|24680+130357.2|25119|157129.1|24992+157129.2|25431|130349.1|24969+130349.2|25408|108931.1|24813+108931.2|25252|10463.1|24976+10463.2|25414|351919.1|24867+351919.2|25305|211660.1|24904+211660.2|25342|14023.1|24989+14023.2|25427|42972.1|24950+42972.2|25387|179435.1|24766+179435.2|25203|171591.1|24765+171591.2|25202|377298.1|24822+377298.2|25259|355804.1|25003+355804.2|25440|63386.1|24927+63386.2|25364|158644.1|24796+158644.2|25233|366714.1|24686+366714.2|25122|70692.1|24884+70692.2|25320|263616.1|24696+263616.2|25132|321384.1|24763+321384.2|25199|209806.1|24768+209806.2|25204|90573.1|24755+90573.2|25191|217829.1|24859+217829.2|25295|216101.1|24985+216101.2|25421|94990.1|25008+94990.2|25444|351730.1|24764+351730.2|25200|33742.1|24923+33742.2|25358|119102.1|24883+119102.2|25318|349708.1|24949+349708.2|25384|182997.1|24817+182997.2|25252|1836.1|24797+1836.2|25232|95961.1|24946+95961.2|25381|287949.1|24764+287949.2|25199|200274.1|24808+200274.2|25243|162078.1|24858+162078.2|25292|175327.1|24898+175327.2|25332|236395.1|25009+236395.2|25443|166630.1|24867+166630.2|25301|355747.1|24880+355747.2|25314|33305.1|24906+33305.2|25340|118738.1|24915+118738.2|25349|99942.1|24757+99942.2|25190|60389.1|24874+60389.2|25307|51695.1|24818+51695.2|25251|82626.1|24869+82626.2|25302|7070.1|24715+7070.2|25147|311865.1|25003+311865.2|25435|98688.1|24921+98688.2|25353|48233.1|24879+48233.2|25310|227372.1|25002+227372.2|25433|267091.1|24815+267091.2|25246|85943.1|24705+85943.2|25136|60398.1|24917+60398.2|25348|288856.1|24981+288856.2|25411|214735.1|24716+214735.2|25146|350983.1|24830+350983.2|25260|236965.1|24706+236965.2|25136|282383.1|24681+282383.2|25111|108410.1|24889+108410.2|25319|148861.1|24999+148861.2|25429|144911.1|24728+144911.2|25157|180584.1|24845+180584.2|25274|29892.1|24774+29892.2|25202|85584.1|24781+85584.2|25209|43298.1|24859+43298.2|25287|210066.1|24775+210066.2|25203|115368.1|24702+115368.2|25130|335097.1|24711+335097.2|25139|126284.1|24952+126284.2|25380|338201.1|24710+338201.2|25137|187656.1|24835+187656.2|25261|16101.1|24702+16101.2|25128|13531.1|25000+13531.2|25426|160675.1|24739+160675.2|25165|230563.1|24716+230563.2|25141|279844.1|24868+279844.2|25293|234844.1|24762+234844.2|25187|345844.1|24924+345844.2|25349|164640.1|24931+164640.2|25356|250469.1|24973+250469.2|25397|29906.1|24760+29906.2|25184|286995.1|24758+286995.2|25182|119600.1|25001+119600.2|25425|8039.1|24727+8039.2|25151|139677.1|24957+139677.2|25381|222544.1|24882+222544.2|25306|196014.1|24749+196014.2|25172|286975.1|24969+286975.2|25392|325356.1|24840+325356.2|25263|127126.1|24979+127126.2|25402|73682.1|24698+73682.2|25121|186799.1|24978+186799.2|25401|189766.1|24987+189766.2|25409|258326.1|24953+258326.2|25374|285241.1|24982+285241.2|25403|219761.1|24979+219761.2|25399|290253.1|24980+290253.2|25400|215512.1|24708+215512.2|25128|228455.1|24919+228455.2|25338|131492.1|24868+131492.2|25287|134144.1|24717+134144.2|25135|254119.1|24898+254119.2|25315|384191.1|24727+384191.2|25143|115855.1|24801+115855.2|25214|87399.1|24703+87399.2|25116|198184.1|24950+198184.2|25361|366258.1|24782+366258.2|25190|234080.1|25020+234080.2|25426|91207.1|24812+91207.2|25218|235276.1|24924+235276.2|25327";
	return list.find(name) != string::npos;
}

bool containedInContig(const FastQSequence &seq, const vector<ContigInfo> &contigs)
{
	int pos = getReadPosition(seq.Comment);
	if (pos < 0)
		return true;
	int n = contigs.size();
	int readLen = seq.Nucleotides.length();
	//bool show = inList(seq.Comment);
	for (int i = 0; i < n; i++)
	{
		if (contigs[i].Position < 0)
			continue;
		//if (show)
		//{
			//printf("[%i %i) must be in [%i %i)\n", pos, pos + readLen, contigs[i].Position, contigs[i].Position + contigs[i].Length);
			//if (pos >= positions[i] && pos <= positions[i] + len)
			//	return true;
			//if (pos + readLen >= positions[i] && pos + readLen <= positions[i] + len)
			//	return true;
		//}
		if (pos >= contigs[i].Position && pos + readLen <= contigs[i].Position + contigs[i].Length)
			return true;
	}
	return false;
}

bool processPairedReads(const PairedInput &input, const vector<ContigInfo> &contigs)
{
	bool success = true;
	FastQReader left, right;
	FastQWriter leftOut, rightOut;
	if (!left.Open(input.LeftFileName) || !right.Open(input.RightFileName))
		success = false;
	if (success && (!leftOut.Open(input.OutputPrefix + "_1.fastq") || !rightOut.Open(input.OutputPrefix + "_2.fastq")))
		success = false;
	if (success)
	{
		FastQSequence leftSeq, rightSeq;
		bool leftRead = left.Read(leftSeq), rightRead = right.Read(rightSeq);
		while (leftRead && rightRead)
		{	
			//if (inList(leftSeq.Comment) || inList(rightSeq.Comment))
			//	cout << "Left:" << endl;
			//bool leftC = containedInContig(leftSeq, contigs);
			//if (inList(leftSeq.Comment) || inList(rightSeq.Comment))
			//	cout << (leftC ? "in" : "out") << endl;
			//if (inList(leftSeq.Comment) || inList(rightSeq.Comment))
			//	cout << "Right:" << endl;
			//bool rightC = containedInContig(rightSeq, contigs);
			//if (inList(leftSeq.Comment) || inList(rightSeq.Comment))
			//	cout << (rightC ? "in" : "out") << endl;
			if (containedInContig(leftSeq, contigs) && containedInContig(rightSeq, contigs))
				leftOut.Write(leftSeq), rightOut.Write(rightSeq);
			leftRead = left.Read(leftSeq), rightRead = right.Read(rightSeq);
		}
		if (leftRead || rightRead)
			success = false;
	}

	left.Close();
	right.Close();
	leftOut.Close();
	rightOut.Close();
	return success;
}

bool processPairedReads(const vector<PairedInput> &input, const vector<ContigInfo> &contigs)
{
	int n = input.size();
	for (int i = 0; i < n; i++)
	{
		cerr << "    [i] Processing paired input: " << input[i].OutputPrefix << endl;
		if (!processPairedReads(input[i], contigs))
			return false;
	}
	return true;
}

int main(int argc, char *argv[])
{
	if (config.ProcessCommandLine(argc, argv))
	{
		if (!readContigs(config.InputFileName, contigs))
		{
			cerr << "[-] Unable to read contigs: " << config.InputFileName << endl;
			return -2;
		}
		cerr << "[+] Successfully read contigs." << endl;
		cerr << "[i] Processing paired reads:" << endl;
		if (!processPairedReads(config.PairedFilter, contigs))
		{
			cerr << "[-] Unable to process paired reads." << endl;
			return -3;
		}
		cerr << "[+] Successfully processed paired reads." << endl;
		return 0;
	}
	cerr << config.LastError;
	return -1;
}
