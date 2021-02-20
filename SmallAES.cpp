#include "SmallAES.h"
using namespace std;

SmallAES::SmallAES() {
	for (int i = 0; i < 16; i++) {
		state[i / 4][i % 4] = 0;
	}
	mulTable.clear();
	mulTableInver.clear();
	mulTable.resize(mulTableSIZE);
	mulTableInver.resize(mulTableSIZE);
	computeMulTable();
}

UINT32 SmallAES::multiplicationX(UINT32 nibble) {
	UINT32 bitTemp;

	bitTemp = (nibble >> 3) & 0x1;
	nibble = (nibble << 1) & 0xf;

	if (bitTemp == 0)
		return nibble;
	else
		return (nibble ^ 0x3);
}

void SmallAES::computeMulTable() {
	//the input is (x0,x1,x2,x3) in GF(2^16) M*(x0,x1,x2,x3)^T=(y0,y1,y2,y3)^T
	UINT32 x[4],y[4],sum;
	for (UINT32 t = 0; t < 0x10000; t++) {
		x[3] = t & 0xf;
		x[2] = (t >> 4) & 0xf;
		x[1] = (t >> 8) & 0xf;
		x[0] = (t >> 12) & 0xf;

		y[0] = multiplicationX(x[0]) ^ multiplicationX(x[1]) ^ x[1] ^ x[2] ^ x[3];
		y[1] = x[0] ^ multiplicationX(x[1]) ^ multiplicationX(x[2]) ^ x[2] ^ x[3];
		y[2] = x[0] ^ x[1] ^ multiplicationX(x[2]) ^ multiplicationX(x[3]) ^ x[3];
		y[3] = multiplicationX(x[0]) ^ x[0] ^ x[1] ^ x[2] ^ multiplicationX(x[3]);

		sum = y[3] + y[2] * 0x10 + y[1] * 0x100 + y[0] * 0x1000;
		mulTable[t] = sum;
		mulTableInver[sum] = t;
	}
}

UINT32 SmallAES::getMul(UINT32 input) {
	return mulTable[input];
}

UINT32 SmallAES::getMulInver(UINT32 input) {
	return mulTableInver[input];
}

void SmallAES::roundFun(UINT32 input[][4], UINT32 output[][4]) {
	UINT32 dig[4];
	for (int i = 0; i < 4; i++) {
		dig[i] = 0;
		dig[i] = S_SB[input[0][i]]*0x1000 + S_SB[input[1][(1 + i) % 4]]*0x100 + S_SB[input[2][(2 + i) % 4]]*0x10 + S_SB[input[3][(3 + i) % 4]];//SB->SR
		dig[i] = getMul(dig[i]);//matrixMul
		//set dig[i] as col[i] of output
		output[0][i] = (dig[i] >> 12) & 0xf;
		output[1][i] = (dig[i] >> 8) & 0xf;
		output[2][i] = (dig[i] >> 4) & 0xf;
		output[3][i] = dig[i] & 0xf;
	}
}

void SmallAES::roundFunInver(UINT32 input[][4], UINT32 output[][4]) {
	UINT32 dig[4];
	for (int i = 0; i < 4; i++) {
		dig[i] = 0;
		dig[i] = input[0][i] * 0x1000 + input[1][i] * 0x100 + input[2][i] * 0x10 + input[3][i];//SB->SR
		dig[i] = getMulInver(dig[i]);//matrixMul
		
		output[0][i] = (dig[i] >> 12) & 0xf;
		output[1][(i+1)%4] = (dig[i] >> 8) & 0xf;
		output[2][(i+2)%4] = (dig[i] >> 4) & 0xf;
		output[3][(i+3)%4] = dig[i] & 0xf;
	}
	for (int i = 0; i < 16; i++) {
		output[i / 4][i % 4] = S_SBInv[output[i / 4][i % 4]];
	}
}
