#ifndef EDR_READER_H
#define EDR_READER_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "FileSerializer.h"

/* This number should be increased whenever the file format changes! */
static constexpr int edr_version = 5;

struct EdrData
{
    EdrData() = default;

    std::string         name; // group name
    std::string         unit; // unit, such as bar, KJ/mol, etc.
    std::vector<double> value;
};

enum EnumEnx
{
    enxOR,     // Time and ensemble averaged data for orientation restraints
    enxORI,    // Instantaneous data for orientation restraints
    enxORT,    // Order tensor(s) for orientation restraints
    enxDISRE,  // Distance restraint blocks
    enxDHCOLL, // Data about the free energy blocks in this frame
    enxDHHIST, // BAR histogram
    enxDH,     // BAR raw delta H data
    enxNR      // Total number of extra blocks in the current code,
               // note that the enxio code can read files written by
               // future code which contain more blocks.
};

enum DataType
{
    xdr_int,
    xdr_float,
    xdr_double,
    xdr_int64,
    xdr_char,
    xdr_string
};

struct SubBlock
{
    int      nr        = 0;
    DataType type      = xdr_float;
    int      val_alloc = 0;
};

struct Block
{
    EnumEnx               id   = enxOR;
    int                   nsub = 0;
    std::vector<SubBlock> sub;
    int                   nsub_alloc = 0;

    void add_subblocks(int n)
    {
        nsub = n;
        if (n > nsub_alloc)
        {
            for (int i = 0; i < n - nsub_alloc; i++)
            {
                sub.emplace_back(SubBlock());
            }
            nsub_alloc = n;
        }
    }
};

struct Frame
{
    double             t            = 0; // frame time in ps
    double             dt           = 0; // frame dt in ps
    int64_t            step         = 0;
    int64_t            nsteps       = 0;
    int                nsum         = 0;
    int                nre          = 0;
    int                nblock       = 0;
    int                nblock_alloc = 0;
    int                e_size       = 0;
    int                e_alloc      = 0;
    std::vector<Block> block;

    void clear()
    {
        block.clear();
        t = dt = 0;
        step = nsteps = nsum = nre = nblock = nblock_alloc = e_size = e_alloc = 0;
    }

    void add_blocks(int n)
    {
        nblock = n;
        if (n > nblock_alloc)
        {
            for (int i = 0; i < n - nblock_alloc; i++)
            {
                block.emplace_back(Block());
            }
            nblock_alloc = n;
        }
    }
};

/* \brief A class to read gromacs edr binary file */
class EdrReader
{
public:
    EdrReader(const std::string& fname)
        : edr_(std::make_unique<FileSerializer>(fname.c_str(), "rb"))
    {
        do_parser();
        // write_data("ene.csv");
    }
    ~EdrReader() { msg("Finished!\n"); }

    // get all energys
    std::map<std::string, std::vector<double>> get_ene() const;

    //! write data to csv file
    void write_data(const std::string& fout) const;

private:
    //! read edr version magic and nre & call do_edr_strings
    bool do_enexnms();
    //! do group names and units
    bool do_edr_strings();
    //! do parser
    bool do_parser();
    //! do a frame data
    bool do_enx();
    //! do header
    bool do_eheader(int nre_test);

private:
    std::vector<EdrData>            data_;
    std::vector<double>             times_; // all time in ps
    std::unique_ptr<FileSerializer> edr_;
    int                             file_version_ = 0;
    int                             nre_          = 0;     //< the number of items
    bool                            is_old_       = false; // is old format edr
    int                             precision_    = 4;     // precision for read float
    Frame                           fr_;                   // a frame data structure
};

#endif // !EDR_READER_H
