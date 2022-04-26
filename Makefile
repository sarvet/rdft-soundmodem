all:     bin2sym307pm8a mod-pm8a Wyman1x-demod-decode flt-clip-scale-wav

CFLAGS=  -O2 -Wall


bin2sym307pm8a: bin2sym307pm8a.c code1_def.h rs8-4.h rs306-xxx-rgm.h sym_code.h sym_def.h mode_def.h
	gcc $(CFLAGS) -o bin2sym307pm8a -lm bin2sym307pm8a.c


Wyman1x-demod-decode: pm8a-demod-decode6q.o dcom-tf.o dfft3f.o decode2.o
	gcc $(CFLAGS) -o Wyman1x-demod-decode -lm pm8a-demod-decode6q.o dcom-tf.o dfft3f.o decode2.o

pm8a-demod-decode6q.o: pm8a-demod-decode6q.c dcom-t.h exfun-t.h dcom.h sym_code.h pm8a-frame.h mode_def.h pm8a-ref-spec2.h gf307-inva.h
	gcc $(CFLAGS) -c -o pm8a-demod-decode6q.o pm8a-demod-decode6q.c



mod-pm8a: mod-pm8a.o dcom-tf.o
	gcc $(CFLAGS) -o mod-pm8a -lm mod-pm8a.o dcom-tf.o

mod-pm8a.o: mod-pm8a.c dcom-t.h exfun-t.h dcom.h sym_def.h
	gcc $(CFLAGS) -c -o mod-pm8a.o mod-pm8a.c



flt-clip-scale-wav: flt-clip-scale-wav.c
	gcc $(CFLAGS) -o flt-clip-scale-wav -lm flt-clip-scale-wav.c



dcom-tf.o: dcom-tf.c dcom-t.h dcom.h
	gcc $(CFLAGS) -c -o dcom-tf.o dcom-tf.c

dfft3f.o: dfft3f.c 
	gcc $(CFLAGS) -c -o dfft3f.o dfft3f.c

decode2.o: decode2.c dcom.h code2-we.h code2-vectors-bias.h code2.h mode-id-vectors-bias.h rs8-4.h sym_def.h
	gcc $(CFLAGS) -c -o decode2.o decode2.c

clean:
	rm -f *.o bin2sym307pm8a mod-pm8a Wyman1x-demod-decode flt-clip-scale-wav


dcom-t.h: dcom.h
sym_code.h: sym_def.h
sym_def.h: code1_def.h
	