#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>     

int lat,lng,alt;

unsigned int get_gps_info(int sec);

void main(int argc,char *argv[])
{
	FILE *in, *out;
	char *line_buf = NULL;
	size_t line_buf_size = 0;
	int line_count = 0;
	ssize_t line_size;
	char line[80];
	char file_name[80];
	char letter;
	char array[80];
	int array_num;
	char name[80];
	int year,month,mday,hour,minute,second,ns;
	int start_minute,start_second;
	int old_second = -1;
	int old_minute = -1;
	int i;
	char *pch;
	int md,ab;
	int micro,nano;   // Micro counter, Nano counter
	unsigned short status[9];
	unsigned short data[3];
	int version = 13;
	unsigned short gps_info;
	unsigned char id;
	unsigned short trig_count=0;
    struct tm mytime = {0,0,0,0,0,0,0,0,0};
	time_t rawtime,datatime,endtime;
    struct tm *ptm;

	/* Work in UTC */
    setenv("TZ","UTC",1);
	tzset();
	if (argc != 2) {
		fprintf(stderr,"Usage:  %s filename\n",argv[0]);
		exit(0);
	}
	in = fopen(argv[1],"r");
	if (!in) {
		fprintf(stderr, "Error opening file '%s' for reading\n", argv[1]);
	}
	line_size = getline(&line_buf, &line_buf_size, in);
	letter = line_buf[0];
	id = letter - 'A' + 1;
	line_size = getline(&line_buf, &line_buf_size, in);
	strcpy(array,line_buf);
	array[line_size-1] = '\0';
	line_size = getline(&line_buf, &line_buf_size, in);
	strcpy(name,line_buf);
	name[line_size-1] = '\0';
	line_size = getline(&line_buf, &line_buf_size, in);
	array_num = atoi(line_buf);
	line_size = getline(&line_buf, &line_buf_size, in);
	lat = atof(line_buf)*324000000.0/90.0;
	line_size = getline(&line_buf, &line_buf_size, in);
	lng = atof(line_buf)*648000000.0/180.0;
	line_size = getline(&line_buf, &line_buf_size, in);
	alt = atof(line_buf)*100.0;
	
	/* Get first peak time */
	line_size = getline(&line_buf, &line_buf_size, in);
	year = atoi(line_buf+2);
	month = atoi(line_buf+5);
	mday = atoi(line_buf+8);
	hour = atoi(line_buf+11);
	minute = atoi(line_buf+14);
	second = atoi(line_buf+17);
	mytime.tm_year = year+100;
	mytime.tm_mon = month;
	mytime.tm_mday = mday;
	mytime.tm_hour = hour;
	mytime.tm_min = minute;
	mytime.tm_sec = second;
	datatime = mktime(&mytime); /* Time of first peak in unix seconds */

	/* Start data file at ten minute interval preceeding first peak */
	mytime.tm_min = start_minute = (minute/10)*10;
	mytime.tm_sec = start_second = 0;
	rawtime = mktime(&mytime);
	ptm = gmtime(&rawtime);
	file_name[0] = 'L';
	file_name[1] = letter;
	file_name[2] = '_';
	file_name[3] = '\0';
	strcat(file_name,array);
	strcat(file_name,"_");
	strcat(file_name,name);
	strcat(file_name,"_");
	sprintf(line,"%02d%02d%02d_%02d%02d%02d.dat",ptm->tm_year-100,ptm->tm_mon,ptm->tm_mday,ptm->tm_hour,ptm->tm_min,ptm->tm_sec);
	strcat(file_name,line);
	printf("%s\n",file_name);
	out = fopen(file_name,"w");

	rawtime = rawtime - 1;
	ptm = gmtime(&rawtime);
    gps_info = get_gps_info(second - 1);
    status[0] = 0x8000 | (version << 7) | (ptm->tm_year - 100);
    status[1] = 0x8000 | ((gps_info >> 2) & 0x2000) | 0x20;
    status[2] = 0x8000 | (ptm->tm_sec << 6) | ptm->tm_min;
    status[3] = 0xc000 | (ptm->tm_hour<<9) | (ptm->tm_mday << 4) | ptm->tm_mon;
    status[4] = 0x8000 | trig_count;
    status[5] = 0x8000 | (id <<8) | array_num;
    status[6] = 0x8000;
    status[7] = 0x8000 | (gps_info&0x7fff);
    status[8] = 0xaa55;
    fwrite(status,sizeof(short),9,out);
	rawtime = rawtime + 1;
	endtime = rawtime + 600;


	trig_count = 0;
	for (;;) {
		while (rawtime < datatime) {
			ptm = gmtime(&rawtime);
		    gps_info = get_gps_info(second - 1);
		    status[0] = 0x8000 | (version << 7) | (ptm->tm_year - 100);
		    status[1] = 0x8000 | ((gps_info >> 2) & 0x2000) | 0x20;
		    status[2] = 0x8000 | (ptm->tm_sec << 6) | ptm->tm_min;
		    status[3] = 0xc000 | (ptm->tm_hour<<9) | (ptm->tm_mday << 4) | ptm->tm_mon;
		    status[4] = 0x8000 | trig_count;
		    status[5] = 0x8000 | (id <<8);
		    status[6] = 0x8000;
		    status[7] = 0x8000 | (gps_info&0x7fff);
		    status[8] = 0xaa55;
		    fwrite(status,sizeof(short),9,out);
			trig_count = 0;
			rawtime = rawtime + 1;
		}
		year = atoi(line_buf+2);
		month = atoi(line_buf+5);
		mday = atoi(line_buf+8);
		hour = atoi(line_buf+11);
		minute = atoi(line_buf+14);
		second = atoi(line_buf+17);
		mytime.tm_year = year+100;
		mytime.tm_mon = month;
		mytime.tm_mday = mday;
		mytime.tm_hour = hour;
		mytime.tm_min = minute;
		mytime.tm_sec = second;
		datatime = mktime(&mytime); /* Time of first peak in unix seconds */
		if (datatime > rawtime) continue;

		ns = atoi(line_buf+20);
		micro = ns/10000;
		nano = (ns - micro*10000)/40;
		pch = strtok(line_buf," ");
		pch = strtok(NULL," ");
		md = atoi(pch);
		pch = strtok(NULL," ");
		ab = atoi(pch);
		if (ab > 199) ab = 199;
		data[0] = nano | ((ab&0x0f)<<11) | ((micro & 0x1c000) >> 6);
		data[1] = 0xc000 | (micro & 0x3fff);
		data[2] = md | ((ab&0x3f0)<<4);
		fwrite(data,sizeof(short),3,out);
		trig_count++;

		line_size = getline(&line_buf, &line_buf_size, in);
		if (line_size == -1) break;
	}
	while (rawtime < endtime) {
		ptm = gmtime(&rawtime);
	    gps_info = get_gps_info(second - 1);
	    status[0] = 0x8000 | (version << 7) | (ptm->tm_year - 100);
	    status[1] = 0x8000 | ((gps_info >> 2) & 0x2000) | 0x20;
	    status[2] = 0x8000 | (ptm->tm_sec << 6) | ptm->tm_min;
	    status[3] = 0xc000 | (ptm->tm_hour<<9) | (ptm->tm_mday << 4) | ptm->tm_mon;
	    status[4] = 0x8000 | trig_count;
	    status[5] = 0x8000 | (id <<8);
	    status[6] = 0x8000;
	    status[7] = 0x8000 | (gps_info&0x7fff);
	    status[8] = 0xaa55;
	    fwrite(status,sizeof(short),9,out);
		rawtime = rawtime + 1;
		trig_count = 0;

	}

	fclose(in);
	fclose(out);
}

unsigned int get_gps_info(int sec) {
	unsigned short gps_info;

	switch (sec % 12) {
		case 0:
			gps_info = (lat >> 16) & 0xffff;
			break;
		case 1:
			gps_info = lat & 0xffff;
			break;
		case 2:
			gps_info = (lng >> 16) & 0xffff;
			break;
		case 3:
			gps_info = lng & 0xffff;
			break;
		case 4:
			gps_info = (alt >> 16) & 0xffff;
			break;
		case 5:
			gps_info = alt & 0xffff;
			break;
		case 6:
			gps_info = 0;
			break;
		case 7:
			gps_info = 0;
			break;
		case 8:
			gps_info = 0;
			break;
		case 9:
			gps_info = 0x0707;
			break;
		case 10:
			gps_info = 0x8b09;
			break;
		case 11:
			gps_info = 0x4e4d;
			break;
	}
	return(gps_info);

}
