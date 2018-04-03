
#include <src/projector.h>
#include <src/backprojector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/ctf.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/ml_model.h>
#include <src/exp_model.h>
#include <src/healpix_sampling.h>

class compare_fsc
{
public:
    FileName fn_map1, fn_map2, fn_out;
    IOParser parser;

    void usage()
    {
        parser.writeUsage(std::cerr);
    }

    void read(int argc, char **argv)
    {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("Options");
        fn_map1 = parser.getOption("--map1", "The first input map to calculate FSC");
        fn_map2 = parser.getOption("--map2", "The second input map to calculate FSC");
        fn_out = parser.getOption("--o", "Output FSC XML file name");

        if (parser.checkForErrors())
            REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
        
    }

    void compare()
    {
        Image<DOUBLE> vol1, vol2;
        MultidimArray<DOUBLE> FSC;

        std::cerr << " Reading map1: " << fn_map1 << std::endl;
        vol1.read(fn_map1);
        std::cerr << " Done reading map1! " << std::endl;

        std::cerr << " Reading map2: " << fn_map2 << std::endl;
        vol2.read(fn_map2);
        std::cerr << " Done reading map2! " << std::endl;

        getFSC(vol1(), vol2(), FSC);
        std::cout << " The FSC curve is: " << FSC << std::endl;
        
    }

};


int main(int argc, char *argv[])
{
    time_config();
    compare_fsc cpf;

    try
    {
        cpf.read(argc, argv);

        cpf.compare();
    }

    catch (RelionError XE)
    {
        cpf.usage();
        std::cout << XE;
        exit(1);
    }

    return 0;
}