/*
 * nodeprobe.c -- 'ping' the mond at each argv.  Print results.
 *
 * Improvement:
 *  Since we send the nodename in the initial 'ping', it's easy to keep track
 *  of which nodes have responded.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/time.h>
#include <signal.h>

#include "mond.h"
int verbose = 1;

int sockfd;
char buf[MAXBUFFER];

/* The child process calls do_send */
void do_send(int numnodes, char **nodenames){
  struct sockaddr_in their_addr; /* connector's address information */
  struct hostent *he;
  size_t sin_size;
  int i;

  /* Send out a request to everybody */
  for (i=0; i<numnodes; i++) {
    int len;
    /* NEEDED!  A better way of getting the node names! */
    errno = 0;
    if ((he=gethostbyname(nodenames[i])) == NULL) {  /* get the host info */
      perror("gethostbyname");
      continue;
    }

    memset(&their_addr, 0, sizeof(their_addr)); /* start with zeros */
    their_addr.sin_family = AF_INET;         /* host byte order */
    their_addr.sin_port = htons(PORT);     /* short, network byte order */
    their_addr.sin_addr = *((struct in_addr *)he->h_addr);

    sin_size = sizeof(their_addr);
    errno = 0;
    len = strlen(nodenames[i]) + 1;
    if( len > MAXBUFFER ){
      fprintf(stderr, "name too long: %s\n", nodenames[i]);
      continue;
    }
    if( sendto(sockfd, nodenames[i], len, 0, &their_addr, sizeof(their_addr)) != len ){
      perror( "sendto" );
      continue;
    }
  }
}

char **nodenames_;
int numnodes_;
void
onalarm(int sig){
  int i;
  for(i=0; i<numnodes_; i++){
    if( strlen(nodenames_[i]) > 0 ){
      printf("%s -1\n", nodenames_[i]);
    }
  }
  exit(0);
}

/* The parent process calls do_recv */
void
do_recv(int numnodes, char **nodenames){
  int nleft;
  /* If the alarm goes off, the default action is to terminate the
     process and write 'Alarm clock' on stderr.  That's probably fine.

     There's always the possibility of the alarm clock going off
     during a 'write', though.  Presumably it's possible to block
     signals while writing, but we are are writing very short buffers
     anyway, so it's unlikely (impossible?) that we would get a
     partial write.  In any event, whatever parses the output should
     probably be defensive about partial or garbled data anyway.

     An alarm handler is also necessary if we want to try to figure
     out who hasn't responded.  But that could equally well be handled
     by the 'parser' */
  numnodes_ = numnodes;
  nodenames_ = nodenames;
  signal(SIGALRM, onalarm);
  alarm(5);
  nleft = numnodes;
  while( nleft > 0 ){
    int len;
    int numbytes;
    int j;
    char *ptr;

    errno = 0;
    if ((numbytes=recvfrom(sockfd, buf, MAXBUFFER, 0, NULL, NULL)) < 0) {
      /* On wealtheow I get lots of errno=ECONNREFUSED here.  Why?
	 This errno is not listed in the man page for recvfrom.
	 Nevertheless, it doesn't seem to be fatal to subsequent processing,
	 so exiting is too drastic.  Let's try continue... */
      perror("recv");
      continue;			/* was exit(4)? */
    }

    /* The last byte should be null, but just in case... */
    buf[numbytes-1] = '\0';
    len = strlen(buf);
    errno = 0;
    if( write(1, buf, len) != len ){
      /* Write failures are probably persistent, and hence unrecoverable. */
      perror("write");
      exit(5);
    }

    /* record which one we heard from... */
    /* look for the first space in buf. */
    ptr = strchr(buf, ' ');
    *ptr = '\0';
    for(j=0; j<numnodes; j++){
      if( strcmp(buf, nodenames[j]) == 0 ){
	nodenames[j] = "";
	break;
      }
    }

    nleft--;
  }
  exit(0);
}

int main(int argc, char *argv[])
{
  int numnodes;
  pid_t pid;

  if (argc < 2) {
    fprintf(stderr,"usage: client node [morenodes]\n");
    exit(1);
  }
  /* Read node names from arg list.  It would probably be better to
     get it from a file.  stdin perhaps? */
  numnodes = argc-1;

  errno = 0;
  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
    perror("socket");
    exit(1);
  }
  /* This works even though we don't call bind().  How? */

  errno = 0;
  pid = fork();
  if( pid < 0 ){
    /* error */
    perror("fork");
    exit(2);
  }else if( pid == 0 ){
    /* child */
    do_send(numnodes, argv+1);
  }else{
    /* parent */
    do_recv(numnodes, argv+1);
  }
  exit(0);
}

