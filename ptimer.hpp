#ifndef _xns_ptimer_
#define _xns_ptimer_

#include <cassert>
#include <chrono>
#include <cstdio>
#include <string>
#include "mpi.h"

namespace xns{

/* Brief
 *  This object is a parallel timer designed for distributed-memory performance analysis. It does not deal yet
 *  with OpenMP multithreaded codes. Please, use it carefully to prevent data races or blocks.
 *
 * Definition
 *  In scientific computing, the performance of the codes is a crucial characteristic. The aim of this object
 *  is to provide developers with a complete but rather simple tool to evaluate the performance of the codes
 *  at user-specified control points. Every control point is transformed into a channel within ptimer, whose
 *  hierarchy within the code is also stored.
 *
 *  By design, the only dependencies of ptimer are MPI and OpenMP, two programming standards for shared- and
 *  distributed-memory parallelism.
 *
 *  Two methods are provided in ptimer class to control channels in user-specified control points: start and
 *  pause.
 *   · start(name): given a channel name, the method searches its channel list to find whether the chanel
 *    exists, or not. Thus, inserting a control point for the first time will create a new channel, while
 *    reaching the same control point multiple times will just increase the channel's calls counter.
 *
 *   · pause(name): given a channel name, the method searches its channel list to find whether the channel
 *    exists and is open, or not. If the channel exists and is open, the time, flops and memory traffic that
 *    ocurred after starting the channel will be accumulated. Otherwise, the application will crash.
 *
 *  Inserting a control point within an existing channel will create a new channel and consider it a nested
 *  child of the existing channel. Besides, if the name given to a nested control point is also used outside
 *  its parent, it will create a new channel with a different hierarchy. For instance, following the code
 *  below, the output will be given as shown in the right:
 *
 *    1. ptimer::start(one point)        |   PTIMER
 *    2.   *some calculations*           |    one point........INFO
 *    3. ptimer::start(another)          |     another.........INFO
 *    4.   *more calculations*           |    another..........INFO
 *    5. ptimer::pause(another)          |
 *    6. ptimer::pause(one point)        |
 *    7. ptimer::start(another)          |
 *    8.   *more stuff*                  |
 *    9. ptimer::pause(another)          |
 *
 *  This allows for fine-grained evaluation of kernels and functions. For instance, inside SpMV kernel the
 *  user may be interested in breaking up its total elapsed time in update and compute times.
 *
 *  Channels accumulate the timing information by default. To store additional information such as the number
 *  of operations or the memory traffic, the information must be added manually inside the control points. To
 *  do so, two functions are provided: add_flops and and_bytes, which increment the variables flops and bytes,
 *  respectively. When calling start, the channel will store the current value of these variables to, later
 *  in pause, calculate their increment. For instance, the following code will account for 10 flops and 80
 *  bytes in the control point:
 *
 *    1. ptimer::start(mypoint)
 *    2.   ptimer::add_flops(10)
 *    3.   ptimer::add_bytes(80)
 *    4. ptimer::pause(mypoint)
 *
 *  Users are encouraged to use the methods above to analyze their kernels by counting or estimating the
 *  performance and throughput within a control point.
 *
 * Pending
 *  Extend the parallel timer to allow for multithreaded performance analysis. Calculate the overhead of
 *  storing one channel per thread.
 *
 *  Extend the member method ptimer::print to write in output files a detailed report process by process. A
 *  collective CSV file would be great.
 */
class ptimer{
public:
    ptimer();
    ~ptimer();

    void start(const char* const name);
    void pause(const char* const name);
    void print();
    void reset();

    template<class T>
    void add_flops(const T _flops){
        flops += (double)_flops;
    }
    template<class T>
    void add_bytes(const T _bytes){
        bytes += (double)_bytes;
    }

private:
    ptimer(const ptimer& rhs) = delete;
    ptimer& operator=(const ptimer& rhs) = delete;

     /* This struct is designed to store all the required information of a channel (times, operation count,
      * memory traffic) as an array within ptimer.
     */
    struct chan_t{
        chan_t() : timer(0.0), flops(0.0), bytes(0.0), level(0), outer(-1), calls(0), open(0), name(0){
            name = new char[MAX_CHANNEL_LEN];
        }
        chan_t(const chan_t& rhs) = delete;
        chan_t& operator=(const chan_t& rhs) = delete;
        ~chan_t(){
            if(name){
                delete[] name;
                name = 0;
            }
        }

        double timer; // time variables
        double flops; // flop variables
        double bytes; // data variables
        int level; // channel's nested level
        int outer; // id of parent channel
        int calls; // number of calls to channel
        bool open; // status
        char* name; // channel's name
    };

    static const int MAX_CHANNEL_NUM = 256; // default max number of channels
    static const int MAX_CHANNEL_LEN = 32; // default max length of channel name

    double flops; // floating-point operations counter
    double bytes; // bytes counter
    int nch; // number of channels created
    int nop; // current number of open channels
    int openid; // current open channel id
    chan_t* channels; // array of channels

    double gettime();
    int create(const char* const name);
    int search(const char* const name);
    bool compare(const char* const name, const int ch);
};

/* Default constructor.
 */
ptimer::ptimer() : flops(0.0), bytes(0.0), nch(0), nop(0), openid(-1), channels(0){
    channels = new chan_t[MAX_CHANNEL_NUM];
}

/* Because dynamic memory is used, the destructor of ptimer must ensure the deletion of all channels and its
 * internal members.
 */
ptimer::~ptimer(){
    if(channels){
        delete[] channels;
        channels = 0;
    }
}

/* The member method ptimer::start resumes an existing channel or creates a new one. To determine whether a
 * chanel exists, both the given name and the current open channel, openid, are evaluated. This allows to use
 * the same name in different nested regions. If the channel name given is the same as the currently open
 * channel, ptimer will crash.
 */
void ptimer::start(const char* const name){
    /* to prevent nesting channels with the same name */
    if(openid>=0){
        assert(std::strcmp(channels[openid].name, name) && "channel is already open");
    }
    int ch = this->search(name);

    /* if channel is not found, create new channel */
    if(ch<0){
        ch = this->create(name);
    }

    nop++;
    openid = ch; // last opened channel

    /* initial time is evaluated at the end to prevent overheads from the search routine */
    channels[ch].open = true;
    channels[ch].flops -= flops;
    channels[ch].bytes -= bytes;
    channels[ch].timer -= this->gettime();
}

/* The member method ptimer::pause pauses an existing channel. To determine whether a chanel exists, both the
 * given name and the current open channel, openid, are evaluated. In the case a channel is not found, or is
 * found but was already closed, ptimer will crash.
 */
void ptimer::pause(const char* const name){
    assert(openid>=0 && "no cannel open");
    assert(!std::strcmp(channels[openid].name, name) && "channel is not open");

    /* store the time, flop and data variables */
    channels[openid].timer += this->gettime();
    channels[openid].bytes += bytes;
    channels[openid].flops += flops;
    channels[openid].open = false;
    channels[openid].calls += 1;

    /* if channel is not found, ptimer::pause will crash */
    openid = channels[openid].outer; // the outer of the current channel becomes the current open channel
    nop--;
}

/* The member method ptimer::print prints a report of the timing on the terminal. The output is adjusted to 80
 * characters width.
 */
void ptimer::print(){
    assert(nch && "no channels found");

    double global = 0.0;

    std::vector<double> tmax(nch), tmin(nch), tavg(nch), favg(nch), bavg(nch);

    for(int i=0; i<nch; ++i){
        tmax[i] = channels[i].timer;
        tmin[i] = channels[i].timer;
        tavg[i] = channels[i].timer;
        favg[i] = channels[i].flops;
        bavg[i] = channels[i].bytes;
    }

    MPI_Allreduce(MPI_IN_PLACE, tmax.data(), nch, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, tmin.data(), nch, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, tavg.data(), nch, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, favg.data(), nch, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, bavg.data(), nch, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int world = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &world);

    for(int i=0; i<nch; ++i){
        tavg[i] /= world;
        favg[i] /= world;
        bavg[i] /= world;

        if(channels[i].level==0){
            global += tavg[i];
        }
    }

    int rank = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(!rank){
        printf("--------------------------------------------------------------------------------\n");
        printf(" PARALLEL TIMER REPORT\n");
        printf("--------------------------------------------------------------------------------\n");
        printf(" Channel                             N GFLOPS   GB/s   tavg   tmax   tmin      %%\n");
        printf("--------------------------------------------------------------------------------\n");

        for(int i=0; i<nch; ++i){
            std::string pname = channels[i].name;

            for(int j=0; j<channels[i].level; ++j){
                pname = " " + pname;
            }

            printf("%-31s %6d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                    pname.c_str(),
                    channels[i].calls,
                    favg[i]/tavg[i]/1e9,
                    bavg[i]/tavg[i]/1e9,
                    tavg[i],
                    tmax[i],
                    tmin[i],
                    tavg[i]/(global>1e-30 ? global : 1.0)*100);
        }
    }
}

/* To reset all channels.
 */
void ptimer::reset(){
    if(channels){
        delete[] channels;
        channels = 0;
    }
    channels = new chan_t[MAX_CHANNEL_NUM];
}

/* The std::chrono is a portable solution for time evaluations. Since we are only interested in elapsed time,
 * the method below provides a floating point value with respect to epoch. We accumulate differences of such a
 * value.
 */
double ptimer::gettime(){
    auto now = std::chrono::system_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::duration<double>>(now).count();
}

/* To create new channels.
 */
int ptimer::create(const char* const name){
    assert(name && "channel name pointer is null");
    assert(name[0] && "channel name is empty");
    assert(std::strlen(name)<MAX_CHANNEL_LEN && "channel name is too long");
    assert(nch<MAX_CHANNEL_NUM && "reached maximum number of channels");

    std::strcpy(channels[nch].name, name);

    /* store the already open channel into outer list; if no channel is open, outer is -1 */
    channels[nch].level = nop;
    channels[nch].outer = openid;

    return(nch++);
}

/* To search a given channel.
 */
int ptimer::search(const char* const name){
    for(int ch=(openid + 1); ch<nch; ch++){
        if(this->compare(name, ch)){
            return ch;
        }
    }

    return -1;
}

/* To compare a given channel.
 */
bool ptimer::compare(const char* const name, const int ch){
    return (std::strcmp(channels[ch].name, name)==0 && channels[ch].outer==openid);
}

}
#endif
