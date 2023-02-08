#include "mechanism.H"
const int rmap[724] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723};

// Returns 0-based map of reaction order
void GET_RMAP
(int * _rmap)
{
for (int j=0; j<724; ++j) {
_rmap[j] = rmap[j];
}
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int * i, int * nspec, int * ki, int * nu)
{
const int ns[724] =
     {3,3,2,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,3,3,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,3,3,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,3,3,5,5,4,4,4,4,4,4,2,2,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,3,3,5,4,4,4,4,4,4,4,4,3,4,4,4,4,4,3,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,2,2,4,4,4,4,5,5,4,4,4,4,4,4,4,4,5,5,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,3,3,3,3,5,4,4,4,4,4,5,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,5,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5};
const int kiv[3620] =
     {8,0,9,0,0,14,0,15,0,0,8,13,0,0,0,0,3,10,0,0,4,11,0,0,0,17,0,14,0,0,16,0,17,0,0,6,2,7,0,0,55,12,0,0,0,14,16,1,0,0,56,8,6,0,0,54,6,0,0,0,6,0,54,0,0,5,0,4,0,0,0,4,5,0,0,12,0,54,0,0,0,54,12,0,0,4,0,2,0,0,0,2,4,0,0,3,2,0,0,0,2,3,0,0,0,1,0,0,0,0,0,1,0,0,0,12,6,1,0,0,6,1,12,0,0,60,12,8,0,0,12,8,60,0,0,31,8,3,0,0,8,3,31,0,0,60,20,0,0,0,20,0,60,0,0,21,38,0,0,0,38,0,21,0,0,9,0,8,1,0,8,1,9,0,0,9,4,8,5,0,8,5,9,4,0,9,2,8,4,0,8,4,9,2,0,13,8,15,9,0,15,9,13,8,0,54,4,6,5,0,6,5,54,4,0,6,4,7,0,0,7,0,6,4,0,0,3,2,4,0,2,4,0,3,0,1,2,0,4,0,0,4,1,2,0,5,2,4,0,0,4,5,2,0,0,1,4,0,5,0,0,5,1,4,0,11,4,5,10,0,5,10,11,4,0,14,2,8,54,0,8,54,14,2,0,13,0,15,1,0,15,1,13,0,0,15,3,14,10,0,14,10,15,3,0,13,4,15,5,0,15,5,13,4,0,13,2,15,4,0,15,4,13,2,0,8,10,55,4,0,55,4,8,10,0,6,10,7,4,0,7,4,6,10,0,6,3,7,2,0,7,2,6,3,0,0,54,6,1,0,6,1,0,54,0,54,2,6,4,0,6,4,54,2,0,12,4,5,54,0,5,54,12,4,0,12,0,1,54,0,1,54,12,0,0,12,2,54,4,0,54,4,12,2,0,8,4,12,1,0,12,1,8,4,0,8,2,12,0,0,12,0,8,2,0,8,3,55,2,0,55,2,8,3,0,12,8,9,54,0,9,54,12,8,0,8,54,9,6,0,9,6,8,54,0,10,2,3,4,0,3,4,10,2,0,54,10,12,3,0,12,3,54,10,0,55,3,12,10,0,12,10,55,3,0,8,10,9,3,0,9,3,8,10,0,54,3,6,10,0,6,10,54,3,0,0,10,4,0,0,4,0,10,0,0,0,10,1,3,0,1,3,0,10,0,10,4,5,3,0,5,3,10,4,0,11,3,10,0,0,10,11,3,0,0,0,11,5,4,0,5,4,0,11,0,9,10,8,11,0,8,11,9,10,0,12,10,11,54,0,11,54,12,10,0,17,15,14,0,0,14,17,15,0,0,14,0,17,1,0,17,1,14,0,0,14,4,17,5,0,17,5,14,4,0,16,3,19,4,0,19,4,16,3,0,17,3,16,10,0,16,10,17,3,0,11,2,10,4,0,10,4,11,2,0,16,2,0,19,0,0,19,16,2,0,16,4,18,0,0,18,0,16,4,0,18,0,8,6,0,8,6,18,0,0,18,2,19,4,0,19,4,18,2,0,18,4,5,19,0,5,19,18,4,0,18,0,1,19,0,1,19,18,0,0,19,4,54,0,0,0,19,65,6,0,65,6,0,19,0,19,2,6,0,0,13,3,15,10,0,15,10,13,3,0,13,10,15,11,0,15,11,13,10,0,17,8,16,9,0,16,9,17,8,0,15,8,14,9,0,14,9,15,8,0,17,0,16,1,0,16,1,17,0,0,15,0,8,0,0,8,15,0,0,0,17,3,12,54,0,12,54,17,3,0,13,15,0,0,0,15,0,13,0,0,14,8,17,9,0,17,9,14,8,0,20,8,54,0,0,8,54,20,0,0,20,3,56,10,0,56,10,20,3,0,20,4,56,5,0,56,5,20,4,0,20,0,56,1,0,56,1,20,0,0,20,2,56,4,0,56,4,20,2,0,20,10,56,11,0,56,11,20,10,0,8,20,56,9,0,56,9,8,20,0,37,16,8,0,0,16,8,37,0,0,22,17,8,0,0,17,8,22,0,0,16,8,21,0,0,21,0,16,8,0,22,37,0,0,0,37,0,22,0,0,22,2,18,8,0,18,8,0,22,2,22,2,15,54,0,15,54,22,2,0,22,10,37,11,0,37,11,22,10,0,22,4,37,5,0,37,5,22,4,0,23,17,0,0,0,17,23,0,0,0,23,4,15,18,0,15,18,23,4,0,23,4,37,12,0,37,12,23,4,0,23,4,17,20,0,17,20,23,4,0,23,2,14,18,0,14,18,23,2,0,23,2,21,12,0,21,12,23,2,0,14,3,17,10,0,17,10,14,3,0,24,14,8,0,0,14,8,24,0,0,24,22,0,0,0,22,0,24,0,0,24,3,22,10,0,22,10,24,3,0,22,2,37,4,0,37,4,22,2,0,22,0,37,1,0,37,1,22,0,0,22,0,14,8,0,14,8,22,0,0,29,14,37,0,0,29,23,8,0,0,23,8,29,0,0,25,23,0,0,0,23,0,25,0,0,25,17,14,0,0,17,14,25,0,0,23,26,25,0,0,25,23,26,0,0,25,8,23,9,0,23,9,25,8,0,37,25,22,23,0,22,23,37,25,0,25,3,23,10,0,23,10,25,3,0,25,0,23,1,0,23,1,25,0,0,15,25,13,23,0,13,23,15,25,0,15,25,14,26,0,14,26,15,25,0,17,25,14,23,0,14,23,17,25,0,26,37,8,0,0,37,8,26,0,0,26,17,15,0,0,17,15,26,0,0,26,25,0,0,0,25,0,26,0,0,26,8,25,9,0,25,9,26,8,0,26,0,25,1,0,25,1,26,0,0,26,4,25,5,0,25,5,26,4,0,37,26,22,25,0,22,25,37,26,0,26,10,25,11,0,25,11,26,10,0,26,3,25,10,0,25,10,26,3,0,58,14,15,0,0,14,15,58,0,0,58,26,0,0,0,26,0,58,0,0,3,58,26,10,0,26,10,3,58,0,27,18,8,0,0,18,8,27,0,0,59,15,6,0,0,15,6,59,0,0,28,0,59,1,0,59,1,28,0,0,28,2,59,4,0,59,4,28,2,0,28,4,59,5,0,59,5,28,4,0,28,8,59,9,0,59,9,28,8,0,28,10,59,11,0,59,11,28,10,0,15,28,59,13,0,59,13,15,28,0,28,15,54,0,0,15,54,28,0,0,28,3,59,10,0,59,10,28,3,0,17,28,14,59,0,14,59,17,28,0,28,37,59,22,0,59,22,28,37,0,30,15,37,0,0,15,37,30,0,0,30,0,29,1,0,29,1,30,0,0,30,2,29,4,0,29,4,30,2,0,30,2,54,58,0,30,2,56,24,0,30,4,29,5,0,29,5,30,4,0,30,4,12,58,0,30,4,20,24,0,30,8,29,9,0,29,9,30,8,0,0,11,1,10,0,1,10,0,11,0,54,2,7,0,0,7,0,54,2,0,56,0,18,1,0,18,1,56,0,0,56,2,18,4,0,18,4,56,2,0,8,56,18,9,0,18,9,8,56,0,14,2,57,0,0,57,0,14,2,0,15,2,20,0,0,20,0,15,2,0,17,14,23,0,0,23,0,17,14,0,60,3,20,10,0,20,10,60,3,0,11,3,10,0,0,10,11,3,0,0,17,3,57,2,0,57,2,17,3,0,32,55,4,0,0,55,4,32,0,0,64,3,6,0,19,12,31,32,54,0,32,54,12,31,0,14,31,17,32,0,17,32,14,31,0,31,9,8,32,0,8,32,31,9,0,15,10,60,4,0,60,4,15,10,0,8,31,55,0,0,15,31,60,55,0,31,10,32,3,0,32,3,31,10,0,11,4,5,10,0,5,10,11,4,0,31,55,3,0,0,13,31,15,32,0,15,32,13,31,0,20,31,56,32,0,56,32,20,31,0,33,17,6,0,0,17,6,33,0,0,34,4,33,5,0,33,5,34,4,0,34,0,33,1,0,33,1,34,0,0,34,2,33,4,0,33,4,34,2,0,34,10,33,11,0,33,11,34,10,0,34,8,33,9,0,33,9,34,8,0,34,31,33,32,0,33,32,34,31,0,61,34,0,0,0,34,0,61,0,0,61,17,12,0,0,17,12,61,0,0,61,3,34,10,0,34,10,61,3,0,37,10,61,4,0,61,4,37,10,0,37,31,61,55,0,22,31,37,32,0,37,32,22,31,0,8,4,65,5,0,65,5,8,4,0,62,17,20,0,0,17,20,62,0,0,62,34,8,0,0,34,8,62,0,0,25,10,62,4,0,62,4,25,10,0,25,31,62,55,0,26,4,12,24,0,12,24,26,4,0,26,2,22,12,0,22,12,26,2,0,26,2,14,20,0,14,20,26,2,0,26,2,15,56,0,15,56,26,2,0,26,4,15,20,0,15,20,26,4,0,26,4,13,56,0,13,56,26,4,0,26,2,59,8,0,59,8,26,2,0,26,4,28,8,0,28,8,26,4,0,26,4,59,9,0,59,9,26,4,0,26,31,25,32,0,25,32,26,31,0,35,63,3,0,0,63,3,35,0,0,35,39,4,0,0,39,4,35,0,0,63,22,12,4,0,36,63,0,0,0,63,36,0,0,0,28,31,59,32,0,59,32,28,31,0,28,25,59,26,0,59,26,28,25,0,21,10,14,6,4,14,6,4,21,10,21,10,38,11,0,38,11,21,10,0,22,3,37,10,0,37,10,22,3,0,22,8,37,9,0,37,9,22,8,0,15,22,13,37,0,13,37,15,22,0,37,10,17,12,4,17,12,4,37,10,37,0,21,1,0,21,1,37,0,0,37,8,21,9,0,21,9,37,8,0,15,37,13,21,0,13,21,15,37,0,15,37,14,22,0,14,22,15,37,0,17,37,14,21,0,14,21,17,37,0,21,22,37,0,0,37,21,22,0,0,37,3,34,4,0,34,4,37,3,0,21,3,38,10,0,38,10,21,3,0,38,0,64,1,0,64,1,38,0,0,21,4,38,5,0,38,5,21,4,0,21,2,14,6,0,14,6,21,2,0,64,4,16,54,0,16,54,64,4,0,37,21,0,0,0,21,0,37,0,0,21,0,38,1,0,38,1,21,0,0,21,8,38,9,0,38,9,21,8,0,21,37,38,22,0,38,22,21,37,0,38,4,64,5,0,64,5,38,4,0,38,3,18,54,0,18,54,38,3,0,36,3,58,0,0,3,58,36,0,0,57,18,0,0,0,18,0,57,0,0,57,3,12,6,4,39,57,20,4,0,37,3,21,10,0,21,10,37,3,0,37,3,57,12,0,57,12,37,3,0,37,3,16,12,4,19,3,7,54,0,7,54,19,3,0,8,3,12,4,0,12,4,8,3,0,14,1,8,0,0,8,14,1,0,0,40,3,10,66,0,10,66,40,3,0,40,4,5,66,0,5,66,40,4,0,0,40,1,66,0,1,66,0,40,0,40,2,66,4,0,66,4,40,2,0,10,40,11,66,0,11,66,10,40,0,8,40,9,66,0,9,66,8,40,0,31,40,32,66,0,32,66,31,40,0,66,6,24,0,0,6,24,66,0,0,41,15,18,0,0,15,18,41,0,0,42,18,24,0,0,18,24,42,0,0,43,3,10,44,0,10,44,43,3,0,43,4,5,44,0,5,44,43,4,0,0,43,1,44,0,1,44,0,43,0,43,2,44,4,0,44,4,43,2,0,10,43,11,44,0,11,44,10,43,0,8,43,9,44,0,9,44,8,43,0,31,43,32,44,0,32,44,31,43,0,44,6,58,0,0,6,58,44,0,0,65,9,8,0,0,8,65,9,0,0,13,65,15,8,0,15,8,13,65,0,65,3,6,0,4,65,1,8,0,0,8,0,65,1,0,65,2,6,0,0,65,4,12,0,0,12,0,65,4,0,65,7,12,6,0,12,6,65,7,0,65,8,14,0,0,14,0,65,8,0,65,18,14,6,0,14,6,65,18,0,45,67,0,0,0,67,0,45,0,0,45,68,0,0,0,68,0,45,0,0,45,69,0,0,0,69,0,45,0,0,45,70,0,0,0,70,0,45,0,0,45,24,58,0,0,24,58,45,0,0,0,45,67,1,0,67,1,0,45,0,0,45,68,1,0,68,1,0,45,0,0,45,69,1,0,69,1,0,45,0,0,45,70,1,0,70,1,0,45,0,45,2,67,4,0,67,4,45,2,0,45,2,68,4,0,68,4,45,2,0,45,2,69,4,0,69,4,45,2,0,45,2,70,4,0,70,4,45,2,0,45,4,67,5,0,67,5,45,4,0,45,4,68,5,0,68,5,45,4,0,45,4,69,5,0,69,5,45,4,0,45,4,70,5,0,70,5,45,4,0,10,45,67,11,0,67,11,10,45,0,10,45,68,11,0,68,11,10,45,0,10,45,69,11,0,69,11,10,45,0,10,45,70,11,0,70,11,10,45,0,8,45,67,9,0,67,9,8,45,0,8,45,68,9,0,68,9,8,45,0,8,45,69,9,0,69,9,8,45,0,8,45,70,9,0,70,9,8,45,0,45,3,67,10,0,67,10,45,3,0,45,3,68,10,0,68,10,45,3,0,45,3,69,10,0,69,10,45,3,0,45,3,70,10,0,70,10,45,3,0,15,45,13,67,0,13,67,15,45,0,15,45,13,68,0,13,68,15,45,0,15,45,13,69,0,13,69,15,45,0,15,45,13,70,0,13,70,15,45,0,17,45,14,67,0,14,67,17,45,0,17,45,14,68,0,14,68,17,45,0,17,45,14,69,0,14,69,17,45,0,17,45,14,70,0,14,70,17,45,0,31,45,67,32,0,67,32,31,45,0,31,45,68,32,0,68,32,31,45,0,31,45,69,32,0,69,32,31,45,0,31,45,70,32,0,70,32,31,45,0,67,45,68,45,0,68,45,67,45,0,67,45,69,45,0,69,45,67,45,0,67,45,70,45,0,70,45,67,45,0,68,45,69,45,0,69,45,68,45,0,68,45,70,45,0,70,45,68,45,0,69,45,70,45,0,70,45,69,45,0,68,22,58,0,0,22,58,68,0,0,68,71,0,0,0,71,0,68,0,0,69,26,24,0,0,26,24,69,0,0,69,71,0,0,0,71,0,69,0,0,69,72,0,0,0,72,0,69,0,0,70,15,30,0,0,15,30,70,0,0,70,72,0,0,0,72,0,70,0,0,68,3,71,10,0,71,10,68,3,0,69,3,71,10,0,71,10,69,3,0,69,3,72,10,0,72,10,69,3,0,70,3,72,10,0,72,10,70,3,0,67,69,0,0,0,69,67,0,0,0,67,70,0,0,0,70,67,0,0,0,68,69,0,0,0,69,68,0,0,0,67,68,0,0,0,68,67,0,0,0,71,4,28,58,0,72,4,28,58,0,71,2,30,20,0,72,2,30,20,0,71,25,24,0,0,25,24,71,0,0,72,25,24,0,0,25,24,72,0,0,46,67,3,0,0,67,3,46,0,0,47,68,3,0,0,68,3,47,0,0,48,69,3,0,0,69,3,48,0,0,49,70,3,0,0,70,3,49,0,0,46,73,0,0,0,73,46,0,0,0,47,74,0,0,0,74,47,0,0,0,47,75,0,0,0,75,47,0,0,0,48,76,0,0,0,76,48,0,0,0,48,77,0,0,0,77,48,0,0,0,48,78,0,0,0,78,48,0,0,0,49,79,0,0,0,79,49,0,0,0,49,50,0,0,0,50,49,0,0,0,74,71,10,0,0,71,10,74,0,0,76,71,10,0,0,71,10,76,0,0,77,72,10,0,0,72,10,77,0,0,50,72,10,0,0,72,10,50,0,0,73,51,4,0,0,75,52,4,0,0,79,52,4,0,0,75,30,20,4,0,78,28,26,4,0,79,22,40,4,0,80,73,3,0,0,73,3,80,0,0,81,75,3,0,0,75,3,81,0,0,82,78,3,0,0,78,3,82,0,0,83,79,3,0,0,79,3,83,0,0,80,84,4,0,0,84,4,80,0,0,81,85,4,0,0,85,4,81,0,0,82,86,4,0,0,86,4,82,0,0,83,87,4,0,0,87,4,83,0,0,84,57,43,4,0,85,27,40,4,0,86,28,41,4,0,87,20,42,4,0,52,4,30,56,5,51,4,14,5,44,52,4,22,5,66,52,10,30,56,11,51,10,14,11,44,52,10,22,11,66};
const int nuv[3620] =
     {-1,-1,1,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,1,0,0,-2,1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,2,0,0,0,-2,1,0,0,0,-1,2,0,0,0,-2,1,0,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-2,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,1,0,-1,1,0,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,-1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,1,-1,1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,0,0,-1,-1,1,0,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,1,1,0,-1,1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,1};
if (*i < 1) {
// Return max num species per reaction
*nspec = 5;
} else {
if (*i > 724) {
*nspec = -1;
} else {
*nspec = ns[*i-1];
for (int j=0; j<*nspec; ++j) {
ki[j] = kiv[(*i-1)*5 + j] + 1;
nu[j] = nuv[(*i-1)*5 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real *  P, const amrex::Real *  T, const amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r)
{
int id; // loop counter
amrex::Real c[54]; // temporary storage
amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (id = 0; id < 54; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, *T);

// convert to chemkin units
for (id = 0; id < 724; ++id) {
q_f[id] *= 1.0e-6;
q_r[id] *= 1.0e-6;
}
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };// temperature cache
amrex::Real invT = 1.0 / tc[1];
// compute the Gibbs free energy
amrex::Real g_RT[54];
gibbs(g_RT, tc);
amrex::Real g_RT_qss[34];
gibbs_qss(g_RT_qss, tc);

amrex::Real sc_qss[34];
// Fill sc_qss here
amrex::Real kf_qss[397], qf_qss[397], qr_qss[397];
comp_k_f_qss(tc, invT, kf_qss);
comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);
comp_sc_qss(sc_qss, qf_qss, qr_qss);
comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 1.008000; // H
awt[1] = 15.999000; // O
awt[2] = 12.011000; // C
awt[3] = 14.007000; // N
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int id; // loop counter
int kd = 4; 
// Zero ncf
for (id = 0; id < kd * 54; ++ id) {
 ncf[id] = 0; 
}

// h
ncf[ 0 * kd + 0 ] = 1; // H

// h2
ncf[ 1 * kd + 0 ] = 2; // H

// o
ncf[ 2 * kd + 1 ] = 1; // O

// o2
ncf[ 3 * kd + 1 ] = 2; // O

// oh
ncf[ 4 * kd + 0 ] = 1; // H
ncf[ 4 * kd + 1 ] = 1; // O

// h2o
ncf[ 5 * kd + 0 ] = 2; // H
ncf[ 5 * kd + 1 ] = 1; // O

// co
ncf[ 6 * kd + 2 ] = 1; // C
ncf[ 6 * kd + 1 ] = 1; // O

// co2
ncf[ 7 * kd + 2 ] = 1; // C
ncf[ 7 * kd + 1 ] = 2; // O

// ch3
ncf[ 8 * kd + 2 ] = 1; // C
ncf[ 8 * kd + 0 ] = 3; // H

// ch4
ncf[ 9 * kd + 2 ] = 1; // C
ncf[ 9 * kd + 0 ] = 4; // H

// ho2
ncf[ 10 * kd + 0 ] = 1; // H
ncf[ 10 * kd + 1 ] = 2; // O

// h2o2
ncf[ 11 * kd + 0 ] = 2; // H
ncf[ 11 * kd + 1 ] = 2; // O

// ch2o
ncf[ 12 * kd + 2 ] = 1; // C
ncf[ 12 * kd + 0 ] = 2; // H
ncf[ 12 * kd + 1 ] = 1; // O

// c2h6
ncf[ 13 * kd + 2 ] = 2; // C
ncf[ 13 * kd + 0 ] = 6; // H

// c2h4
ncf[ 14 * kd + 2 ] = 2; // C
ncf[ 14 * kd + 0 ] = 4; // H

// c2h5
ncf[ 15 * kd + 2 ] = 2; // C
ncf[ 15 * kd + 0 ] = 5; // H

// c2h2
ncf[ 16 * kd + 2 ] = 2; // C
ncf[ 16 * kd + 0 ] = 2; // H

// c2h3
ncf[ 17 * kd + 2 ] = 2; // C
ncf[ 17 * kd + 0 ] = 3; // H

// ch2co
ncf[ 18 * kd + 2 ] = 2; // C
ncf[ 18 * kd + 0 ] = 2; // H
ncf[ 18 * kd + 1 ] = 1; // O

// hcco
ncf[ 19 * kd + 2 ] = 2; // C
ncf[ 19 * kd + 0 ] = 1; // H
ncf[ 19 * kd + 1 ] = 1; // O

// ch3cho
ncf[ 20 * kd + 2 ] = 2; // C
ncf[ 20 * kd + 0 ] = 4; // H
ncf[ 20 * kd + 1 ] = 1; // O

// c3h4-a
ncf[ 21 * kd + 2 ] = 3; // C
ncf[ 21 * kd + 0 ] = 4; // H

// c3h6
ncf[ 22 * kd + 2 ] = 3; // C
ncf[ 22 * kd + 0 ] = 6; // H

// c4h6
ncf[ 23 * kd + 2 ] = 4; // C
ncf[ 23 * kd + 0 ] = 6; // H

// nc3h7
ncf[ 24 * kd + 2 ] = 3; // C
ncf[ 24 * kd + 0 ] = 7; // H

// c4h7
ncf[ 25 * kd + 2 ] = 4; // C
ncf[ 25 * kd + 0 ] = 7; // H

// c4h8-1
ncf[ 26 * kd + 2 ] = 4; // C
ncf[ 26 * kd + 0 ] = 8; // H

// ch3coch2
ncf[ 27 * kd + 2 ] = 3; // C
ncf[ 27 * kd + 0 ] = 5; // H
ncf[ 27 * kd + 1 ] = 1; // O

// c2h5cho
ncf[ 28 * kd + 2 ] = 3; // C
ncf[ 28 * kd + 0 ] = 6; // H
ncf[ 28 * kd + 1 ] = 1; // O

// c5h9
ncf[ 29 * kd + 2 ] = 5; // C
ncf[ 29 * kd + 0 ] = 9; // H

// c5h10-1
ncf[ 30 * kd + 2 ] = 5; // C
ncf[ 30 * kd + 0 ] = 10; // H

// ch3o2
ncf[ 31 * kd + 2 ] = 1; // C
ncf[ 31 * kd + 0 ] = 3; // H
ncf[ 31 * kd + 1 ] = 2; // O

// ch3o2h
ncf[ 32 * kd + 2 ] = 1; // C
ncf[ 32 * kd + 0 ] = 4; // H
ncf[ 32 * kd + 1 ] = 2; // O

// c2h3co
ncf[ 33 * kd + 2 ] = 3; // C
ncf[ 33 * kd + 0 ] = 3; // H
ncf[ 33 * kd + 1 ] = 1; // O

// c2h3cho
ncf[ 34 * kd + 2 ] = 3; // C
ncf[ 34 * kd + 0 ] = 4; // H
ncf[ 34 * kd + 1 ] = 1; // O

// c4h8ooh1-3o2
ncf[ 35 * kd + 2 ] = 4; // C
ncf[ 35 * kd + 0 ] = 9; // H
ncf[ 35 * kd + 1 ] = 4; // O

// pc4h9o2
ncf[ 36 * kd + 2 ] = 4; // C
ncf[ 36 * kd + 0 ] = 9; // H
ncf[ 36 * kd + 1 ] = 2; // O

// c3h5-a
ncf[ 37 * kd + 2 ] = 3; // C
ncf[ 37 * kd + 0 ] = 5; // H

// c3h3
ncf[ 38 * kd + 2 ] = 3; // C
ncf[ 38 * kd + 0 ] = 3; // H

// nc4ket13
ncf[ 39 * kd + 2 ] = 4; // C
ncf[ 39 * kd + 0 ] = 8; // H
ncf[ 39 * kd + 1 ] = 3; // O

// nc3h7cho
ncf[ 40 * kd + 2 ] = 4; // C
ncf[ 40 * kd + 0 ] = 8; // H
ncf[ 40 * kd + 1 ] = 1; // O

// c2h5coch2
ncf[ 41 * kd + 2 ] = 4; // C
ncf[ 41 * kd + 0 ] = 7; // H
ncf[ 41 * kd + 1 ] = 1; // O

// nc3h7coch2
ncf[ 42 * kd + 2 ] = 5; // C
ncf[ 42 * kd + 0 ] = 9; // H
ncf[ 42 * kd + 1 ] = 1; // O

// nc4h9cho
ncf[ 43 * kd + 2 ] = 5; // C
ncf[ 43 * kd + 0 ] = 10; // H
ncf[ 43 * kd + 1 ] = 1; // O

// nc4h9co
ncf[ 44 * kd + 2 ] = 5; // C
ncf[ 44 * kd + 0 ] = 9; // H
ncf[ 44 * kd + 1 ] = 1; // O

// nc7h16
ncf[ 45 * kd + 2 ] = 7; // C
ncf[ 45 * kd + 0 ] = 16; // H

// c7h15o2-1
ncf[ 46 * kd + 2 ] = 7; // C
ncf[ 46 * kd + 0 ] = 15; // H
ncf[ 46 * kd + 1 ] = 2; // O

// c7h15o2-2
ncf[ 47 * kd + 2 ] = 7; // C
ncf[ 47 * kd + 0 ] = 15; // H
ncf[ 47 * kd + 1 ] = 2; // O

// c7h15o2-3
ncf[ 48 * kd + 2 ] = 7; // C
ncf[ 48 * kd + 0 ] = 15; // H
ncf[ 48 * kd + 1 ] = 2; // O

// c7h15o2-4
ncf[ 49 * kd + 2 ] = 7; // C
ncf[ 49 * kd + 0 ] = 15; // H
ncf[ 49 * kd + 1 ] = 2; // O

// c7h14ooh4-3
ncf[ 50 * kd + 2 ] = 7; // C
ncf[ 50 * kd + 0 ] = 15; // H
ncf[ 50 * kd + 1 ] = 2; // O

// c7h14o1-3
ncf[ 51 * kd + 2 ] = 7; // C
ncf[ 51 * kd + 0 ] = 14; // H
ncf[ 51 * kd + 1 ] = 1; // O

// c7h14o2-4
ncf[ 52 * kd + 2 ] = 7; // C
ncf[ 52 * kd + 0 ] = 14; // H
ncf[ 52 * kd + 1 ] = 1; // O

// n2
ncf[ 53 * kd + 3 ] = 2; // N

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(4);
ename[0] = "H";
ename[1] = "O";
ename[2] = "C";
ename[3] = "N";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(54);
kname[0] = "h";
kname[1] = "h2";
kname[2] = "o";
kname[3] = "o2";
kname[4] = "oh";
kname[5] = "h2o";
kname[6] = "co";
kname[7] = "co2";
kname[8] = "ch3";
kname[9] = "ch4";
kname[10] = "ho2";
kname[11] = "h2o2";
kname[12] = "ch2o";
kname[13] = "c2h6";
kname[14] = "c2h4";
kname[15] = "c2h5";
kname[16] = "c2h2";
kname[17] = "c2h3";
kname[18] = "ch2co";
kname[19] = "hcco";
kname[20] = "ch3cho";
kname[21] = "c3h4-a";
kname[22] = "c3h6";
kname[23] = "c4h6";
kname[24] = "nc3h7";
kname[25] = "c4h7";
kname[26] = "c4h8-1";
kname[27] = "ch3coch2";
kname[28] = "c2h5cho";
kname[29] = "c5h9";
kname[30] = "c5h10-1";
kname[31] = "ch3o2";
kname[32] = "ch3o2h";
kname[33] = "c2h3co";
kname[34] = "c2h3cho";
kname[35] = "c4h8ooh1-3o2";
kname[36] = "pc4h9o2";
kname[37] = "c3h5-a";
kname[38] = "c3h3";
kname[39] = "nc4ket13";
kname[40] = "nc3h7cho";
kname[41] = "c2h5coch2";
kname[42] = "nc3h7coch2";
kname[43] = "nc4h9cho";
kname[44] = "nc4h9co";
kname[45] = "nc7h16";
kname[46] = "c7h15o2-1";
kname[47] = "c7h15o2-2";
kname[48] = "c7h15o2-3";
kname[49] = "c7h15o2-4";
kname[50] = "c7h14ooh4-3";
kname[51] = "c7h14o1-3";
kname[52] = "c7h14o2-4";
kname[53] = "n2";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(Jac[ 55 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 55 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 55 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 55;
int offset_col = nc * 55;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(Jac[55*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 55*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 55*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
