#include <stdio.h>

#include <fcntl.h>
 #include <unistd.h>
#include <bits/types/struct_timeval.h>

#define SIZE 1024
int main(int argc, char **argv) {
    struct timeval time;
    if(argc!=2){
        perror("wrong param");
        return -3;
    }
    int file_desc;
    if ((file_desc = open(argv[1],O_RDONLY)) == -1){
        perror("can't open file");
        return -2;

    }
    int file_desc1;
    if ((file_desc1 = open("/dev/tty", O_RDONLY | O_NDELAY)) == -1) {
        perror("/dev/tty");
        return -2;
    }

    char c;
    int str_pos = 0;
    int str = 0;
    int str_len[SIZE];
    int offset[SIZE];
    char buf[SIZE];
    while(read(file_desc,&c,1)){/// read byte by byte, symbol == 1 byte
        if(c == '\n') { /// \n ! string!!\n
            ++str_pos;
            str_len[str] = str_pos;
            ++str;
            offset[str] = lseek(file_desc,0,SEEK_CUR);/// offset + cur_pos
            str_pos = 0;

        }
        else ++str_pos;


        }

    for(int i = 0;i < str;++i){
        printf("%d\n",offset[i]);
    }

    int number;
    while(1){

        printf("write number of string");
        scanf("%d",&number);
        if(number == 0){
            break;
        }
        if(number < 0){
            printf("wrong number: ");
            return -1;
        }
        lseek(file_desc,offset[number],SEEK_SET);
        if(read(file_desc,buf,str_len[number])){
            write(STDOUT_FILENO,buf,str_len[number]);
        }else {
            printf("bad number: ");

        }




    }
    close(file_desc);
    close(file_desc1);

    return 0;
    }





