#include "SmallAES.h"
#include <iostream> 
#include <ctime>
using namespace std;

void mulMap(UINT32 m2[], UINT32 m3[]) {
	for (UINT32 i = 0; i < 0x100; i++) {
		UINT32 h = i & 0x80; /* hi bit */
		UINT32 b = (i << 1)&0xff;

		if (h == 0x80) {
			b ^= 0x1b; /* Rijndael's Galois field */
		}
		m2[i] = b;
		m3[i] = b ^ i;
	}
}

void testAES(UINT32 m2[],UINT32 m3[]) {
	int testTimes = 10;
	cout << "total test times:" << testTimes << endl;
	UINT32 table[0x100];
	UINT32 zero[2] = { 0,0 };
	for (int j = 0; j < testTimes; j++) {
		int a0 = rand() & 0xff;
		int a1 = rand() & 0xff;
		int a2 = rand() & 0xff;
		int a3 = rand() & 0xff;
		int a4 = rand() & 0xff;
		int a5 = rand() & 0xff;
		while ((a0 == 0 || a1 == 0 || a2 == 0 || a3 == 0) || (a2==a3)) {
			a0 = rand() & 0xff;
			a1 = rand() & 0xff;
			a2 = rand() & 0xff;
			a3 = rand() & 0xff;
		}
		for (int i = 0; i < 0x100; i++) {
			table[i] = 0;
		}
		cout << "current times:" << j << endl;
		cout << "related constants:" << hex << a0 << " " << a1 << " " << a2 << " " << a3 << " " << a4 << " " << a5 << endl;
		zero[0] = 0;
		zero[1] = 0;
		for (UINT32 i = 0; i <= 0xffffffff;i++) {
			int x0 = i & 0xff;
			int x1 = (i>>8) & 0xff;
			int x2 = (i>>16) & 0xff;
			int x3 = (i>>24) & 0xff;
			//form-1
			int sum0 = 0;
			sum0 = sum0^ m2[SBOX[m2[SBOX[x0]] ^ m3[SBOX[x1]] ^ SBOX[x2] ^ SBOX[x3]]];
			sum0 = sum0^ m2[SBOX[m2[SBOX[x0^a0]] ^ m3[SBOX[x1^a1]] ^ SBOX[x2^a2] ^ SBOX[x3^a3]]];
			sum0 = sum0 ^ x0^a4;
			zero[0] = zero[0] ^ SBOX[sum0];
			//form-2
			int sum1 = 0;
			sum1 = sum1^ SBOX[m2[SBOX[x0]] ^ m3[SBOX[x1]] ^ SBOX[x2] ^ SBOX[x3]];
			sum1 = sum1^ SBOX[m2[SBOX[x0^a0]] ^ m3[SBOX[x1^a1]] ^ SBOX[x2^a2] ^ SBOX[x3^a3]];
			sum1 = sum1 ^ x2 ^ a5;
			zero[1] = zero[1] ^ SBOX[sum1];
			if (i == 0xffffffff) {
				break;
			}
		}
		cout << "form-1 sum: " << zero[0];
		if (zero[0] != 0) {
			cout << " not zero";
		}
		cout << endl;
		cout << "form-2 sum: " << zero[1];
		if (zero[1] != 0) {
			cout << " not zero";
		}
		cout << endl;
	}
}

void testEven() {
	SmallAES s;
	int cnt0[0x10];
	int cnt1[0x10];
	bool isEven[2] = { 1,1 };
	for (int t = 0; t < 0x10000; t++) {
		//cout << "current c1[D(0)]:" << hex << t << endl;
		int t0 = t & 0xf;
		int t1 = (t >> 4) & 0xf;
		int t2 = (t >> 8) & 0xf;
		int t3 = (t >> 12) & 0xf;

		for (int i = 0; i < 0x10; i++) {
			cnt0[i] = 0;
			cnt1[i] = 0;
		}

		for (int a = 0; a < 0x10000; a++) {
			int a0 = a & 0xf;
			int a1 = (a >> 4) & 0xf;
			int a2 = (a >> 8) & 0xf;
			int a3 = (a >> 12) & 0xf;

			int sum0 = s.multiplicationX(S_SB[s.multiplicationX(S_SB[a0]) ^ s.multiplicationX(S_SB[a1]) ^ S_SB[a1] ^ S_SB[a2] ^ S_SB[a3]]);
			sum0 = sum0 ^ s.multiplicationX(S_SB[s.multiplicationX(S_SB[a0 ^ t0]) ^ s.multiplicationX(S_SB[a1 ^ t1]) ^ S_SB[a1 ^ t1] ^ S_SB[a2 ^ t2] ^ S_SB[a3 ^ t3]]);
			sum0 = sum0 ^ a0;
			cnt0[sum0]++;

			int sum1 = S_SB[s.multiplicationX(S_SB[a0]) ^ s.multiplicationX(S_SB[a1]) ^ S_SB[a1] ^ S_SB[a2] ^ S_SB[a3]];
			sum1 = sum1 ^ S_SB[s.multiplicationX(S_SB[a0 ^ t0]) ^ s.multiplicationX(S_SB[a1 ^ t1]) ^ S_SB[a1 ^ t1] ^ S_SB[a2 ^ t2] ^ S_SB[a3 ^ t3]];
			sum1 = sum1 ^ a2;
			cnt1[sum1]++;
		}
		//analyze table
		for (int i = 0; i < 0x10; i++) {
			if (cnt0[i] % 2 != 0) {
				isEven[0] = false;
				cout << 111 << endl;
			}
			if (cnt1[i] % 2 != 0) {
				isEven[1] = false;
				cout << 111 << endl;
			}
		}
	}
	if (isEven[0] && isEven[1]) {
		cout << "The same value always appears even times in all cases" << endl;
	}
	else {
		cout << "the property does not hold for the small-scale aes" << endl;
	}
}

void copyState(UINT32 src[][4], UINT32 des[][4]) {
	for (int t = 0; t < 16; t++) {
		des[t / 4][t % 4] = src[t / 4][t % 4];
	}
}

void addState(UINT32 src[][4], UINT32 des[][4]) {
	for (int t = 0; t < 16; t++) {
		des[t / 4][t % 4] = src[t / 4][t % 4] ^ des[t / 4][t % 4];
	}
}

void clearState(UINT32 state[][4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			state[i][j] = 0;
		}
	}
}

void outputState(UINT32 state[][4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << hex<<state[i][j] << " ";
		}
		cout << endl;
	}
}

void computeWeakKeys(SmallAES &sAES, UINT32 c[][4], int digNum,UINT32 key[][4]) {
	UINT32 dig[4],sum;
	UINT32 can = 0,canre=0;
	UINT32 size = 0x10000;
	bool find[4] = { false,false,false,false};
	bool isweak = true;
	UINT32 input[4][4], output[4][4];

	while (!find[0]) {
		can = rand() % size;
		dig[0] = (can >> 12) & 0xf;
		dig[1] = (can >> 8) & 0xf;
		dig[2] = (can >> 4) & 0xf;
		dig[3] = can & 0xf;

		sum = S_SB[dig[0] ^ c[0][0]] * 0x1000 + S_SB[dig[1] ^ c[1][1]] * 0x100 + S_SB[dig[2] ^ c[2][2]] * 0x10 + S_SB[dig[3] ^ c[3][3]];
		sum = sAES.getMul(sum);
		sum = (sum >> 12) & 0xf;

		if (sum == dig[0]) {
			//a possible weak key
			key[0][0] = dig[0];
			key[1][1] = dig[1];
			key[2][2] = dig[2];
			key[3][3] = dig[3];

			find[0] = true;
		}
	}
	//cout << "dig 0 is fixed" << endl;

	for (int i = 1; i < 4; i++) {
		while (!find[i]) {
			canre = rand() % size;
			dig[0] = (canre >> 12) & 0xf;
			dig[1] = (canre >> 8) & 0xf;
			dig[2] = (canre >> 4) & 0xf;
			dig[3] = canre & 0xf;

			sum = S_SB[dig[0] ^ c[0][i]] * 0x1000 + S_SB[dig[1] ^ c[1][(i + 1) % 4]] * 0x100 + S_SB[dig[2] ^ c[2][(i + 2) % 4]] * 0x10 + S_SB[dig[3] ^ c[3][(i + 3) % 4]];
			sum = sAES.getMul(sum);
			sum = (sum >> 12 - 4 * i) & 0xf;
			if (sum == key[i][i]) {
				key[0][i] = dig[0];
				key[1][(i + 1) % 4] = dig[1];
				key[2][(i + 2) % 4] = dig[2];
				key[3][(i + 3) % 4] = dig[3];
				
				find[i] = true;
			}
		}
		//cout << "dig"<<i <<" is fixed" << endl;
	}
	//check
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			input[i][j] = key[i][j] ^ c[i][j];
		}
	}

	sAES.roundFun(input, output);
	for (int i = 0; i < 4; i++) {
		if (output[i][i] != key[i][i]) {
			isweak = false;
		}
	}
	if (isweak) {
		cout << "a weak key!" << endl;
	}
	else {
		cout << "not a weak key!" << endl;
	}
}

void computeWeakKeysForKeyRecovery(SmallAES& sAES, UINT32 di[], UINT32 c[][4], int digNum, UINT32 key[][4]) {
	UINT32 dig[4], sum;
	UINT32 can = 0, canre = 0;
	UINT32 size = 0x10000;
	bool find[4] = { false,false,false,false };
	bool isweak = true;
	UINT32 input[4][4], output[4][4];

	key[0][0] = di[0];
	key[1][1] = di[1];
	key[2][2] = di[2];
	key[3][3] = di[3];

	//cout << "dig 0 is fixed" << endl;

	for (int i = 1; i < 4; i++) {
		while (!find[i]) {
			canre = rand() % size;
			dig[0] = (canre >> 12) & 0xf;
			dig[1] = (canre >> 8) & 0xf;
			dig[2] = (canre >> 4) & 0xf;
			dig[3] = canre & 0xf;

			sum = S_SB[dig[0]] * 0x1000 + S_SB[dig[1]] * 0x100 + S_SB[dig[2]] * 0x10 + S_SB[dig[3]];
			sum = sAES.getMul(sum);
			sum = (sum >> 12 - 4 * i) & 0xf;
			if (sum == key[i][i]) {
				key[0][i] = dig[0];
				key[1][(i + 1) % 4] = dig[1];
				key[2][(i + 2) % 4] = dig[2];
				key[3][(i + 3) % 4] = dig[3];

				find[i] = true;
			}
		}
		//cout << "dig"<<i <<" is fixed" << endl;
	}
	//check
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			input[i][j] = key[i][j];
		}
	}

	sAES.roundFun(input, output);
	for (int i = 1; i < 4; i++) {
		if (output[i][i] != key[i][i]) {
			isweak = false;
		}
	}
	if (isweak) {
		cout << "a weak key!" << endl;
	}
	else {
		cout << "not a weak key!" << endl;
	}
}

bool isWeakKey(SmallAES& sAES, UINT32 c[][4], int digNum,UINT32 key[][4]) {
	UINT32 input[4][4], output[4][4];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			input[i][j] = key[i][j] ^ c[i][j];
		}
	}

	sAES.roundFun(input, output);
	for (int i = 0; i < 4; i++) {
		if (output[i][i] != key[i][i]) {
			return true;
		}
	}

	return false;
}

void verify4RDistinguisher() {
	SmallAES sAES;
	UINT32 input[4][4],output[4][4],sum[4][4],tmp1[4][4],tmp2[4][4];
	//input[0][0], input[1][1], input[2][2], input[3][3] take all possible values
	for (int i = 0; i < 16; i++) {
		input[i / 4][i % 4] = 0;
		sum[i / 4][i % 4] = 0;
		output[i / 4][i % 4] = 0;
	}

	for (int i = 0; i < 0x10000; i++) {
		input[0][0] = (i >> 12) & 0xf;
		input[1][1] = (i >> 8) & 0xf;
		input[2][2] = (i >> 4) & 0xf;
		input[3][3] = i & 0xf;

		sAES.roundFun(input, tmp1);
		sAES.roundFun(tmp1, tmp2);
		sAES.roundFun(tmp2, tmp1);
		sAES.roundFun(tmp1, output);

		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				sum[r][c] = sum[r][c] ^ output[r][c];
			}
		}
	}

	//output sum
	for (int r = 0; r < 4; r++) {
		for (int c = 0; c < 4; c++) {
			cout << sum[r][c] << " ";
		}
		cout << endl;
	}
}

void AEGISDistinguisher(int isRandomKey) {
	//set seed
	//define constants
	UINT32 c0[4][4], c1[4][4];
	UINT32 iv[4][4], k[4][4];
	//randomly choose a value for c0 and c1
	for (int i = 0; i < 16; i++) {
		c0[i / 4][i % 4] = rand() % 16;
		c1[i / 4][i % 4] = rand() % 16;
	}
	//start the 5-round permutation
	UINT32 s0[4][4], s1[4][4], s2[4][4], s3[4][4], s4[4][4];//s0=k+iv, s1=c0, s2=c2, s3=k+c0, s4=k+c1
	UINT32 ts0[4][4], ts1[4][4], ts2[4][4], ts3[4][4], ts4[4][4];
	UINT32 sum0[4][4], sum1[4][4], sumQ[4][4], sum4[4][4];
	for (int i = 0; i < 16; i++) {
		sum0[i / 4][i % 4] = 0;
		sum1[i / 4][i % 4] = 0;
		sumQ[i / 4][i % 4] = 0;
		sum4[i / 4][i % 4] = 0;
	}

	SmallAES sAES;
	//the 5-round weak-key distinguisher

	//choose a weak key
	if(isRandomKey)
		computeWeakKeys(sAES, c1, 0, k);

	else {
		for (int i = 0; i < 16; i++) {
			k[i / 4][i % 4] = rand() % 0x10;
		}
	}

	if (isWeakKey(sAES, c1, 0, k)) {
		cout << "not a weak key" << endl;
	}

	//vary iv[D(0)]
	for (int i = 0; i < 16; i++) {
		iv[i / 4][i % 4] = rand()%0xf;
	}
	for (UINT32 val = 0; val < 0x10000; val++) {
		
		iv[0][0] = (val >> 12) & 0xf;
		iv[1][1] = (val >> 8) & 0xf;
		iv[2][2] = (val >> 4) & 0xf;
		iv[3][3] = (val) & 0xf;

		//load state
		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				s0[r][c] = iv[r][c] ^ k[r][c];
			}
		}
		copyState(c0, s1);
		copyState(c1, s2);
		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				s3[r][c] = c0[r][c] ^ k[r][c];
				s4[r][c] = c1[r][c] ^ k[r][c];
			}
		}

		//update state
		for (int r = 0; r < 5; r++) {
			sAES.roundFun(s0, ts0);
			sAES.roundFun(s1, ts1);
			sAES.roundFun(s2, ts2);
			sAES.roundFun(s3, ts3);
			sAES.roundFun(s4, ts4);
			if (r % 2 == 0) {
				for (int i = 0; i < 16; i++) {
					s0[i / 4][i % 4] = s0[i / 4][i % 4] ^ k[i / 4][i % 4];
				}
			}
			else {
				for (int i = 0; i < 16; i++) {
					s0[i / 4][i % 4] = s0[i / 4][i % 4] ^ k[i / 4][i % 4] ^ iv[i / 4][i % 4];
				}
			}
			//update state
			addState(ts4, s0);
			addState(ts0, s1);
			addState(ts1, s2);
			addState(ts2, s3);
			addState(ts3, s4);
		}
		/*///after 3-round permutation, check the first column of s1
		cout << s1[0][0] << " " << s1[1][0] << " " << s1[2][0] << " " << s1[3][0] << endl;
		system("pause");
		*/
		///compute the sum of the state after 5-round permutation
		addState(s0, sum0);
		addState(s1, sum1);
		addState(s4, sum4);
		//compute the quadratic part
		for (int i = 0; i < 16; i++) {
			sumQ[i / 4][i % 4] = sumQ[i / 4][i % 4] ^ (s2[i / 4][i % 4] & s3[i / 4][i % 4]);
		}
	}
	//check the sum
	/*cout << "sum 0:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 1; j < 4; j++) {
			cout << sum0[i][j] << " ";
		}
		cout << endl;
	}
	cout << "sum 1:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << sum1[i][j] << " ";
		}
		cout << endl;
	}
	cout << "sum 4:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << sum4[i][j] << " ";
		}
		cout << endl;
	}
	cout << "sum Q:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << sumQ[i][j] << " ";
		}
		cout << endl;
	}*/
	cout << "sum of output:" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << (sumQ[i][j]^sum1[i][j]^sum4[i][j]) << " ";
		}
		cout << endl;
	}
}

void keyRecoveryAEGIS() {
	//define constants
	UINT32 c0[4][4], c1[4][4];
	UINT32 iv[4][4], k[4][4];
	//randomly choose a value for c0 and c1
	for (int i = 0; i < 16; i++) {
		c0[i / 4][i % 4] = rand() % 16;
		c1[i / 4][i % 4] = rand() % 16;
	}
	//start the 5-round permutation
	UINT32 s0[4][4], s1[4][4], s2[4][4], s3[4][4], s4[4][4];//s0=k+iv, s1=c0, s2=c2, s3=k+c0, s4=k+c1
	UINT32 ts0[4][4], ts1[4][4], ts2[4][4], ts3[4][4], ts4[4][4];
	UINT32 sum0[4][4], sum1[4][4], sumQ[4][4], sum4[4][4],zero[4][4];
	for (int i = 0; i < 16; i++) {
		sum0[i / 4][i % 4] = 0;
		sum1[i / 4][i % 4] = 0;
		sumQ[i / 4][i % 4] = 0;
		sum4[i / 4][i % 4] = 0;
		zero[i / 4][i % 4] = 0;
	}

	SmallAES sAES;
	//the 5-round weak-key distinguisher

	//choose a weak key
	computeWeakKeys(sAES, c1, 0, k);

	//not choose a weak key
	/*for (int i = 0; i < 16; i++) {
		k[i / 4][i % 4] = rand() % 0x10;
	}*/

	if (isWeakKey(sAES, c1, 0, k)) {
		cout << "not a weak key" << endl;
	}

	//vary iv[D(0)]
	for (int i = 0; i < 16; i++) {
		iv[i / 4][i % 4] = rand() % 0x10;
	}
	//compute the candidates of k[D(0)] if a weak key is used
	UINT32 dig[4],sum,constantT;
	vector<UINT32> candidates;
	candidates.clear();
	for (int can = 0; can < 0x10000; can++) {
		dig[0] = (can >> 12) & 0xf;
		dig[1] = (can >> 8) & 0xf;
		dig[2] = (can >> 4) & 0xf;
		dig[3] = can & 0xf;

		sum = S_SB[dig[0] ^ c1[0][0]] * 0x1000 + S_SB[dig[1] ^ c1[1][1]] * 0x100 + S_SB[dig[2] ^ c1[2][2]] * 0x10 + S_SB[dig[3] ^ c1[3][3]];
		sum = sAES.getMul(sum);
		sum = (sum >> 12) & 0xf;

		if (sum == dig[0]) {
			//a possible weak key
			candidates.push_back(can);
		}
	}
	cout << "possible number of solutions for k[D(0)]" << ": 0x" << hex << candidates.size() << endl;

	//randomly choose a value for (T[1,2,3][0]) [T=c0 ^ A(k^iv)]
	constantT = rand() % 0x1000;
	UINT32 varyT = 0;
	UINT32 dc = c0[0][0] * 0x1000 + c0[1][0] * 0x100 + c0[2][0] * 0x10 + c0[3][0];
	UINT32 dv = 0;
	UINT32 keyCan[4];

	for (UINT32 it = 0; it < candidates.size(); it++) {
		keyCan[0] = (candidates[it] >> 12) & 0xf;
		keyCan[1] = (candidates[it] >> 8) & 0xf;
		keyCan[2] = (candidates[it] >> 4) & 0xf;
		keyCan[3] = (candidates[it]) & 0xf;
		//correct guess
		/*keyCan[0] = k[0][0];
		keyCan[1] = k[1][1];
		keyCan[2] = k[2][2];
		keyCan[3] = k[3][3];*/

		//initialize sum
		copyState(zero, sum0);
		copyState(zero, sum1);
		copyState(zero, sumQ);
		copyState(zero, sum4);

		for (UINT32 val = 0; val < 16; val++) {
			varyT = constantT + val * 0x1000;
			varyT = varyT ^ dc;
			dv = sAES.getMulInver(varyT);
			iv[0][0] = (dv >> 12) & 0xf;
			iv[1][1] = (dv >> 8) & 0xf;
			iv[2][2] = (dv >> 4) & 0xf;
			iv[3][3] = (dv ) & 0xf;
			iv[0][0] = S_SBInv[iv[0][0]] ^ keyCan[0];
			iv[1][1] = S_SBInv[iv[1][1]] ^ keyCan[1];
			iv[2][2] = S_SBInv[iv[2][2]] ^ keyCan[2];
			iv[3][3] = S_SBInv[iv[3][3]] ^ keyCan[3];

			//load state
			for (int r = 0; r < 4; r++) {
				for (int c = 0; c < 4; c++) {
					s0[r][c] = iv[r][c] ^ k[r][c];
				}
			}
			//cout << "input:" << endl;
			//outputState(iv);

			copyState(c0, s1);
			copyState(c1, s2);
			for (int r = 0; r < 4; r++) {
				for (int c = 0; c < 4; c++) {
					s3[r][c] = c0[r][c] ^ k[r][c];
					s4[r][c] = c1[r][c] ^ k[r][c];
				}
			}

			//update state
			for (int r = 0; r < 5; r++) {
				sAES.roundFun(s0, ts0);
				sAES.roundFun(s1, ts1);
				sAES.roundFun(s2, ts2);
				sAES.roundFun(s3, ts3);
				sAES.roundFun(s4, ts4);
				if (r % 2 == 0) {
					for (int i = 0; i < 16; i++) {
						s0[i / 4][i % 4] = s0[i / 4][i % 4] ^ k[i / 4][i % 4];
					}
				}
				else {
					for (int i = 0; i < 16; i++) {
						s0[i / 4][i % 4] = s0[i / 4][i % 4] ^ k[i / 4][i % 4] ^ iv[i / 4][i % 4];
					}
				}
				//update state
				addState(ts4, s0);
				addState(ts0, s1);
				addState(ts1, s2);
				addState(ts2, s3);
				addState(ts3, s4);
			}
			///compute the sum of the state after 5-round permutation
			addState(s0, sum0);
			addState(s1, sum1);
			addState(s4, sum4);
			//compute the quadratic part
			for (int i = 0; i < 16; i++) {
				sumQ[i / 4][i % 4] = sumQ[i / 4][i % 4] ^ (s2[i / 4][i % 4] & s3[i / 4][i % 4]);
			}
		}
		//check sum
		/*cout << "sum 1:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << sum1[i][j] << " ";
			}
			cout << endl;
		}
		cout << "sum 4:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << sum4[i][j] << " ";
			}
			cout << endl;
		}
		cout << "sum Q:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << sumQ[i][j] << " ";
			}
			cout << endl;
		}*/
		//check whether the guessed key is correct
		bool isKeyCorrect = true;
		for (int i = 0; i < 4; i++) {
			for (int j = 1; j < 4; j++) {
				sum1[i][j] = sum1[i][j] ^ sum4[i][j] ^ sumQ[i][j];
				if (sum1[i][j] != 0) {
					isKeyCorrect = false;
				}
			}
		}
		if (isKeyCorrect) {
			cout << "correct guess!" << endl;
			cout << "correct key: " << hex<<k[0][0] << " " << k[1][1] << " " << k[2][2] << " " << k[3][3] << endl;
			cout << "guessed key:" << hex<<keyCan[0]<< " " << keyCan[1] << " " << keyCan[2] << " " << keyCan[3] << endl;
			//system("pause");
		}
	}
	cout << "all possible values have been traversed" << endl;
}

//analysis of 8-round Tiaoxin
void Tiaoxin8RDistinguisher(int keyChoice) {
	UINT32 z0[4][4], z1[4][4], zero[4][4];//used constants
	UINT32 k[4][4], iv[4][4];
	UINT32 q[4][4], input[4][4], output[4][4];
	//Tiaoxin state
	UINT32 U[3][4][4], W[4][4][4], Y[6][4][4];
	UINT32 TU[3][4][4], TW[4][4][4], TY[6][4][4];
	UINT32 sum[13][4][4],sumQ1[4][4],sumQ2[4][4];
	UINT32 test[4][4];
	for (int i = 0; i < 13; i++) {
		clearState(sum[i]);
	}
	clearState(sumQ1);
	clearState(sumQ2);
	clearState(zero);

	SmallAES sAES;

	for (int i = 0; i < 16; i++) {
		k[i / 4][i % 4] = rand() % 0x10;
		z0[i / 4][i % 4] = rand() % 0x10;
		z1[i / 4][i % 4] = rand() % 0x10;
		input[i / 4][i % 4] = rand() % 0x10;
	}

	if (keyChoice != 1 && keyChoice!=2 && keyChoice!=3) {
		cout << "a random key" << endl;
	}

	//select a weak key (A(k)[D(0)]=z0[D(0)])
	//A(k)[D(0)]=z0[D(0)]
	input[0][0] = z0[0][0];
	input[1][1] = z0[1][1];
	input[2][2] = z0[2][2];
	input[3][3] = z0[3][3];
	if (keyChoice == 1) {
		cout << "a weak key" << endl;
		cout << "Condition: A(k)[D(0)]=z0[D(0)]" << endl;
		sAES.roundFunInver(input, k);
	}

	//A^2(k)[D(0)]=z0[D(0)]
	if (keyChoice == 2) {
		cout << "a weak key" << endl;
		cout << "Condition: A^2(k)[D(0)]=z0[D(0)]" << endl;
		sAES.roundFunInver(input, output);
		sAES.roundFunInver(output, k);
	}

	//A^2(k)[D(0)]=A(k)[D(0)]
	if (keyChoice == 3) {
		cout << "Condition: A^2(k)[D(0)]=A(k)[D(0)]" << endl;
		computeWeakKeys(sAES, zero, 0, output);
		sAES.roundFunInver(output, k);
	}
	//check
	/*sAES.roundFun(k, input);
	sAES.roundFun(input, output);
	for (int i = 0; i < 4; i++) {
		if (input[i][i] != output[i][i]) {
			cout << "wrong" << endl;
		}
	}*/
	
	//iv=A^-1(q)
	for (int i = 0; i < 16; i++) {
		q[i / 4][i % 4] = rand() % 0x10;
	}

	for (UINT32 val = 0; val < 0x10000; val++) {
		q[0][0] = (val >> 12) & 0xf;
		q[1][1] = (val >> 4) & 0xf;
		q[2][2] = (val >> 8) & 0xf;
		q[3][3] = (val) & 0xf;
		//q[0][0] = val;

		sAES.roundFunInver(q, iv);

		//update function of Tiaoxin
		//load state
		copyState(k, U[0]);
		copyState(k, U[1]);
		copyState(iv, U[2]);

		copyState(k, W[0]);
		copyState(k, W[1]);
		copyState(iv, W[2]);
		copyState(z0, W[3]);

		copyState(k, Y[0]);
		copyState(k, Y[1]);
		copyState(iv, Y[2]);
		copyState(z1, Y[3]);
		copyState(zero, Y[4]);
		copyState(zero, Y[5]);

		for (int r = 0; r < 8; r++) {
			//update U
			sAES.roundFun(U[2], TU[2]);
			sAES.roundFun(U[0], TU[0]);
			addState(TU[2], U[0]);
			addState(z0, U[0]);//U[0] = z0 + U[0] + A(U[2])

			copyState(U[1], U[2]);//U[2]=U[1]
			copyState(TU[0], U[1]);//U[1] = A(U[0])
			//update W
			sAES.roundFun(W[3], TW[3]);
			sAES.roundFun(W[0], TW[0]);
			addState(TW[3], W[0]);
			addState(z1, W[0]);
			copyState(W[2], W[3]);
			copyState(W[1], W[2]);
			copyState(TW[0], W[1]);
			//update Y
			sAES.roundFun(Y[5], TY[5]);
			sAES.roundFun(Y[0], TY[0]);
			addState(TY[5], Y[0]);
			addState(z0, Y[0]);
			copyState(Y[4], Y[5]);
			copyState(Y[3], Y[4]);
			copyState(Y[2], Y[3]);
			copyState(Y[1], Y[2]);
			copyState(TY[0], Y[1]);
		}

		//sum up 
		addState(U[0], sum[0]);
		addState(U[1], sum[1]);
		addState(U[2], sum[2]);

		addState(W[0], sum[3]);
		addState(W[1], sum[4]);
		addState(W[2], sum[5]);
		addState(W[3], sum[6]);

		addState(Y[0], sum[7]);
		addState(Y[1], sum[8]);
		addState(Y[2], sum[9]);
		addState(Y[3], sum[10]);
		addState(Y[4], sum[11]);
		addState(Y[5], sum[12]);

		//sum qradratic part
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				sumQ1[i][j] = sumQ1[i][j] ^ (Y[3][i][j] & W[3][i][j]);
				sumQ2[i][j] = sumQ2[i][j] ^ (Y[5][i][j] & U[2][i][j]);
			}
		}
	}
	//inverse the sum[0]
	UINT32 tmp = 0;
	/*for (int i = 0; i < 4; i++) {
		tmp = sum[0][0][i] * 0x1000 + sum[0][1][i] * 0x100 + sum[0][2][i] * 0x10 + sum[0][3][i];
		tmp = sAES.getMulInver(tmp);
		sum[0][0][i] = (tmp >> 12) & 0xf;
		sum[0][1][i] = (tmp >> 8) & 0xf;
		sum[0][2][i] = (tmp >> 4) & 0xf;
		sum[0][3][i] = (tmp) & 0xf;
	}

	for (int i = 0; i < 4; i++) {
		tmp = sum[2][0][i] * 0x1000 + sum[2][1][i] * 0x100 + sum[2][2][i] * 0x10 + sum[2][3][i];
		tmp = sAES.getMulInver(tmp);
		sum[2][0][i] = (tmp >> 12) & 0xf;
		sum[2][1][i] = (tmp >> 8) & 0xf;
		sum[2][2][i] = (tmp >> 4) & 0xf;
		sum[2][3][i] = (tmp) & 0xf;
	}

	cout << "U0:" << endl;
	outputState(sum[0]);
	cout << "U2:" << endl;
	outputState(sum[2]);
	cout << "W1:" << endl;
	outputState(sum[4]);
	cout << "Q1:" << endl;
	outputState(sumQ1);

	cout << "U1:" << endl;
	outputState(sum[1]);*/

	//final sum
	UINT32 finalSum[4][4];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			finalSum[i][j] = sum[0][i][j] ^ sum[2][i][j] ^ sum[4][i][j] ^ sumQ1[i][j];
		}
	}
	//apply the inverse of MC
	for (int i = 0; i < 4; i++) {
		tmp = finalSum[0][i] * 0x1000 + finalSum[1][i] * 0x100 + finalSum[2][i] * 0x10 + finalSum[3][i];
		tmp = sAES.getMulInver(tmp);
		finalSum[0][i] = (tmp >> 12) & 0xf;
		finalSum[1][i] = (tmp >> 8) & 0xf;
		finalSum[2][i] = (tmp >> 4) & 0xf;
		finalSum[3][i] = (tmp) & 0xf;
	}
	//outputState(finalSum);
	bool isOutput[4][4] = {
		0,1,1,1,
		1,1,1,0,
		1,1,0,1,
		1,0,1,1
	};
	cout << "the sum of the targeted bytes:";
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (isOutput[i][j]) {
				cout << finalSum[i][j] << " ";
			}
		}
	}
	cout << endl;
}

//the key-recovery attack on 8-roound Tiaoxin when using a weak constant
void keyRecoveryTiaoxin() {
	UINT32 z0[4][4], z1[4][4], zero[4][4];//used constants
	UINT32 k[4][4], iv[4][4];
	UINT32 q[4][4], input[4][4], output[4][4];
	//Tiaoxin state
	UINT32 U[3][4][4], W[4][4][4], Y[6][4][4];
	UINT32 TU[3][4][4], TW[4][4][4], TY[6][4][4];
	UINT32 sum[13][4][4], sumQ1[4][4], sumQ2[4][4];
	UINT32 test[4][4];
	for (int i = 0; i < 13; i++) {
		clearState(sum[i]);
	}
	clearState(sumQ1);
	clearState(sumQ2);
	clearState(zero);

	SmallAES sAES;

	for (int i = 0; i < 16; i++) {
		k[i / 4][i % 4] = rand() % 0x10;
		z0[i / 4][i % 4] = rand() % 0x10;
		z1[i / 4][i % 4] = rand() % 0x10;
		input[i / 4][i % 4] = rand() % 0x10;
	}

	//condition on z0: A(z0)[0][0] = z0[0][0]
	UINT32 tmpZ0[4][4];
	bool find = false;
	for (int i = 0; i < 4; i++) {
		z1[i][i] = z0[i][i];
	}
	
	//a weak key
	//first let A(K)[D(0)] = z0[D(0)]
	UINT32 di[4];
	di[0] = z0[0][0];
	di[1] = z0[1][1];
	di[2] = z0[2][2];
	di[3] = z0[3][3];
	/*di[0] = rand() % 0x10;
	di[1] = rand() % 0x10;
	di[2] = rand() % 0x10;
	di[3] = rand() % 0x10;*/

	computeWeakKeysForKeyRecovery(sAES, di, zero, 0, output);
	sAES.roundFunInver(output, k);
	for (int i = 0; i < 4; i++) {
		cout << hex<<k[i][i] << " ";
	}
	cout << endl;
	//check the conditions on the key
	sAES.roundFun(k, input);
	for (int i = 0; i < 4; i++) {
		if (input[i][i] != z0[i][i]) {
			cout << "not a weak key" << endl;
		}
	}
	sAES.roundFun(input, output);
	for (int i = 1; i < 4; i++) {
		if (output[i][i] != z0[i][i]) {
			cout << "not a weak key" << endl;
		}
	}

	//key-recovery attack
	//we need guess K[D(0)].
	UINT32 dig[4];
	//set constant part 
	UINT32 GCon[4][4],QCon[4][4];
	UINT32 G[4][4];
	for (int i = 0; i < 16; i++) {
		QCon[i / 4][i % 4] = rand() % 0x10;
		GCon[i / 4][i % 4] = rand() % 0x10;
	}
	//compute the candidates for k[D(0)]
	vector<UINT32> candidates;
	candidates.clear();
	UINT32 tsum = 0;
	for (int can = 0; can < 0x10000; can++) {
		dig[0] = (can >> 12) & 0xf;
		dig[1] = (can >> 8) & 0xf;
		dig[2] = (can >> 4) & 0xf;
		dig[3] = can & 0xf;

		tsum = S_SB[dig[0]] * 0x1000 + S_SB[dig[1]] * 0x100 + S_SB[dig[2]] * 0x10 + S_SB[dig[3]];
		tsum = sAES.getMul(tsum);
		tsum = (tsum >> 12) & 0xf;

		if (tsum == z0[0][0]) {
			//a possible weak key
			candidates.push_back(can);
		}
	}
	cout << "possible number of solutions for k[D(0)]" << ": 0x" << hex << candidates.size() << endl;

	for (UINT32 can = 0; can < candidates.size(); can++) {
		dig[0] = (candidates[can] >> 12) & 0xf;//guessed K[0][0]
		dig[1] = (candidates[can] >> 8) & 0xf;//guessed K[1][1]
		dig[2] = (candidates[can] >> 4) & 0xf;//guessed K[2][2]
		dig[3] = (candidates[can]) & 0xf;//guessed K[3][3]

		//dig[0] = k[0][0];
		//dig[1] = k[1][1];
		//dig[2] = k[2][2];
		//dig[3] = k[3][3];

		//clear sum
		for (int i = 0; i < 13; i++) {
			clearState(sum[i]);
		}
		clearState(sumQ1);
		clearState(sumQ2);

		//generate the inputs
		//G = A(q+K+z0) -> q = A^-1(G) + K + z0
		//set constant for G
		for (int i = 0; i < 16; i++) {
			G[i / 4][i % 4] = GCon[i / 4][i % 4];
		}
		//inverse G
		for (UINT32 it = 0; it < 0x10; it++) {
			G[0][0] = it;
			sAES.roundFunInver(G, q);
			for (int t = 0; t < 4; t++) {
				q[t][t] = q[t][t] ^ dig[t] ^ z0[t][t];
			}
			//randomly choose values for other bytes
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					if (i != j) {
						q[i][j] = QCon[i][j];
					}
				}
			}
			//compute N
			sAES.roundFunInver(q, iv);
			//update function of Tiaoxin
			//load state
			copyState(k, U[0]);
			copyState(k, U[1]);
			copyState(iv, U[2]);

			copyState(k, W[0]);
			copyState(k, W[1]);
			copyState(iv, W[2]);
			copyState(z0, W[3]);

			copyState(k, Y[0]);
			copyState(k, Y[1]);
			copyState(iv, Y[2]);
			copyState(z1, Y[3]);
			copyState(zero, Y[4]);
			copyState(zero, Y[5]);

			for (int r = 0; r < 8; r++) {
				//update U
				sAES.roundFun(U[2], TU[2]);
				sAES.roundFun(U[0], TU[0]);
				addState(TU[2], U[0]);
				addState(z0, U[0]);//U[0] = z0 + U[0] + A(U[2])

				copyState(U[1], U[2]);//U[2]=U[1]
				copyState(TU[0], U[1]);//U[1] = A(U[0])
				//update W
				sAES.roundFun(W[3], TW[3]);
				sAES.roundFun(W[0], TW[0]);
				addState(TW[3], W[0]);
				addState(z1, W[0]);
				copyState(W[2], W[3]);
				copyState(W[1], W[2]);
				copyState(TW[0], W[1]);
				//update Y
				sAES.roundFun(Y[5], TY[5]);
				sAES.roundFun(Y[0], TY[0]);
				addState(TY[5], Y[0]);
				addState(z0, Y[0]);
				copyState(Y[4], Y[5]);
				copyState(Y[3], Y[4]);
				copyState(Y[2], Y[3]);
				copyState(Y[1], Y[2]);
				copyState(TY[0], Y[1]);
			}

			//sum up 
			addState(U[0], sum[0]);
			addState(U[1], sum[1]);
			addState(U[2], sum[2]);

			addState(W[0], sum[3]);
			addState(W[1], sum[4]);
			addState(W[2], sum[5]);
			addState(W[3], sum[6]);

			addState(Y[0], sum[7]);
			addState(Y[1], sum[8]);
			addState(Y[2], sum[9]);
			addState(Y[3], sum[10]);
			addState(Y[4], sum[11]);
			addState(Y[5], sum[12]);

			//sum qradratic part
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					sumQ1[i][j] = sumQ1[i][j] ^ (Y[3][i][j] & W[3][i][j]);
				}
			}
		}

		UINT32 tmp = 0;
		UINT32 finalSum[4][4];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				finalSum[i][j] = sum[0][i][j] ^ sum[2][i][j] ^ sum[4][i][j] ^ sumQ1[i][j];
			}
		}
		//apply the inverse of MC
		for (int i = 0; i < 4; i++) {
			tmp = finalSum[0][i] * 0x1000 + finalSum[1][i] * 0x100 + finalSum[2][i] * 0x10 + finalSum[3][i];
			tmp = sAES.getMulInver(tmp);
			finalSum[0][i] = (tmp >> 12) & 0xf;
			finalSum[1][i] = (tmp >> 8) & 0xf;
			finalSum[2][i] = (tmp >> 4) & 0xf;
			finalSum[3][i] = (tmp) & 0xf;
		}
		//cout << "final sum" <<endl;
		//outputState(finalSum);
		//system("pause");
		bool isOutput[4][4] = {
			0,1,1,1,
			0,1,1,0,
			0,1,0,1,
			0,0,1,1
		};
		bool find = true;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (isOutput[i][j] == 1) {
					if (finalSum[i][j] != 0) {
						find = false;
					}
				}
			}
		}
		if (find) {
			cout << "find the key!" << endl;
			//cout << k[0][0] << " " << k[1][1] << " " << k[2][2] << " " << k[3][3] << endl;
			cout << dig[0] << " " << dig[1] << " " << dig[2] << " " << dig[3] << endl;
			//system("pause");
			//break;
		}
	}
	cout << "all the candidates of k[D(0)] are traversed!" << endl;
}

void realAESTest() {
	UINT32 m2[0x100], m3[0x100];
	mulMap(m2, m3);
	testAES(m2, m3);
}

int main() {
	srand(time(NULL));

	cout << "1 -> Test the wrong integral property on real aes round function (time-consuming)" << endl;
	cout << "2 -> verify the correct integral property on small-scale aes round function" << endl;
	cout << "3 -> 8-round distinguisher for Tiaoxin" << endl;
	cout << "4 -> key recovery for 8-round Tiaoxin with weak constants" << endl;
	cout << "5 -> 5-round distinguisher for AEGIS" << endl;
	cout << "6 -> 5-round key-recovery for AEGIS" << endl;

	int cmd;
	cout << endl << "input command (1/2/3/4/5/6):";
	cin >> cmd;
	if (cmd == 1) {
		realAESTest();
	}
	else if (cmd == 2) {
		testEven();
	}
	else if (cmd == 3) {
		int keyChoice;
		cout << "0 -> a random key" << endl;
		cout << "1/2/3 -> different kinds of weak keys" << endl;
		cout << "please input choice (0/1/2/3): ";
		cin >> keyChoice;
		Tiaoxin8RDistinguisher(keyChoice);
	}
	else if (cmd == 4) {
		for (int i = 0; i < 32; i++) {
			cout << "Times: " << i << endl;
			keyRecoveryTiaoxin();
			cout << endl;
		}
	}
	else if (cmd == 5) {
		int keyChoice;
		cout << "0 -> a random key" << endl;
		cout << "1-> a weak key" << endl;
		cout << "please input choice (0/1): ";
		cin >> keyChoice;
		AEGISDistinguisher(keyChoice);
	}
	else if (cmd == 6) {
		for (int i = 0; i < 32; i++) {
			cout << "Times: " << i << endl;
			keyRecoveryAEGIS();
			cout << endl;
		}
	}
	else {
		cout << "wrong commands" << endl;
	}
	return 0;
}