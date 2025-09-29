#include "EdrReader.h"

#include <fstream>

#include "TprException.h"

#define MAX_STRLEN 1024
constexpr int Vergen = 26; // char use int to read


bool EdrReader::do_enexnms()
{
    int magic;
    if (!edr_->do_int(&magic)) return TPR_FAILED;
    msg("magic= %d\n", magic);
    if (magic > 0)
    {
        // Assume this is an old edr format
        file_version_ = 1;
        nre_          = magic;
        is_old_       = true;
    }
    else
    {
        // new edr format
        is_old_ = false;
        if (magic != -55555)
        {
            THROW_TPR_EXCEPTION(
                "Energy names magic number mismatch, this is not a GROMACS edr file");
        }
        file_version_ = edr_version;
        if (!edr_->do_int(&file_version_)) return TPR_FAILED;
        if (file_version_ > edr_version)
        {
            THROW_TPR_EXCEPTION("Reading tpx file version " + std::to_string(file_version_)
                                + " with version " + std::to_string(edr_version) + " program");
        }
        if (!edr_->do_int(&nre_)) return TPR_FAILED;
    }
    msg("nre= %d\n", nre_);
    msg("file_version_= %d\n", file_version_);

    if (file_version_ != edr_version)
    {
        std::string message = "Note: enx file_version " + std::to_string(file_version_)
                              + ", software version " + std::to_string(edr_version) + "\n";
        msg(message.c_str());
    }

    return do_edr_strings();
}

bool EdrReader::do_edr_strings()
{
    // read items name and unit
    char        name[MAX_STRLEN], unit[MAX_STRLEN];
    std::string gname, gunit;
    data_.resize(nre_);
    for (int i = 0; i < nre_; i++)
    {
        if (!edr_->xdr_string(name, MAX_STRLEN))
        {
            THROW_TPR_EXCEPTION("Error reading edr groups name");
        }

        gname = name;
        if (file_version_ >= 2)
        {
            if (!edr_->xdr_string(unit, MAX_STRLEN))
            {
                THROW_TPR_EXCEPTION("Error reading edr groups unit");
            }
            gunit = unit;
        }
        else { gunit = "kJ/mol"; }
        msg("group %d name= %s, unit= %s\n", i + 1, gname.c_str(), gunit.c_str());

        data_[i].name = gname;
        data_[i].unit = gunit;
    }

    return TPR_SUCCESS;
}

bool EdrReader::do_parser()
{
    int filever = -1;

    // first read
    do_enexnms();

    // loop all frames
    while (do_enx())
    {
        // store fr.t only for fr.nre == nre_
        if (fr_.nre == nre_) { times_.emplace_back(fr_.t); }

        // clear fr data
        fr_.clear();
    }

    return TPR_SUCCESS;
}

bool EdrReader::do_enx()
{
    if (!do_eheader(-1)) return TPR_FAILED;

    // Check sanity of this header
    bool bSane = fr_.nre > 0;
    for (int b = 0; b < fr_.nblock; b++)
    {
        bSane = bSane || (fr_.block[b].nsub > 0);
    }
    if (!(fr_.step >= 0 && bSane))
    {
        THROW_TPR_EXCEPTION("there may be something wrong with energy file");
    }

    if (fr_.nre > fr_.e_alloc)
    {
        for (int i = fr_.e_alloc; i < fr_.nre; i++)
        {
            // TODO:
        }
        fr_.e_alloc = fr_.nre;
    }

    for (int i = 0; i < fr_.nre; i++)
    {
        float ene;
        if (!edr_->do_real(&ene, precision_)) return TPR_FAILED;
        // save energy
        data_[i].value.emplace_back(ene);

        // old edr format
        if (file_version_ == 1 || fr_.nsum > 0 || fr_.nsum > 1)
        {
            float eav, esum, fdum;
            if (!edr_->do_real(&eav, precision_)) return TPR_FAILED;
            if (!edr_->do_real(&esum, precision_)) return TPR_FAILED;

            if (file_version_ == 1)
            {
                // Old, unused real
                if (!edr_->do_real(&fdum, precision_)) return TPR_FAILED;
            }
        }
    }

    // read the blocks
    for (size_t b = 0; b < fr_.nblock; b++)
    {
        for (size_t s = 0; s < fr_.block[b].nsub; s++)
        {
            const int nr = fr_.block[b].sub[s].nr;
            switch (fr_.block[b].sub[s].type)
            {
                case xdr_int:
                {
                    std::vector<int> temp(nr);
                    if (!edr_->do_vector(temp.data(), nr, precision_, Vergen)) return TPR_FAILED;
                    break;
                }
                case xdr_float:
                {
                    std::vector<float> temp(nr);
                    if (!edr_->do_vector(temp.data(), nr, precision_, Vergen)) return TPR_FAILED;
                    break;
                }
                case xdr_double:
                {
                    std::vector<double> temp(nr);
                    if (!edr_->do_vector(temp.data(), nr, precision_, Vergen)) return TPR_FAILED;
                    break;
                }
                case xdr_int64:
                {
                    std::vector<int64_t> temp(nr);
                    if (!edr_->do_vector(temp.data(), nr, precision_, Vergen)) return TPR_FAILED;
                    break;
                }
                case xdr_char:
                {
                    std::vector<unsigned char> temp(nr);
                    if (!edr_->do_vector(temp.data(), nr, precision_, Vergen)) return TPR_FAILED;
                    break;
                }
                case xdr_string:
                {
                    char buff[MAX_STRLEN];
                    for (int i = 0; i < nr; i++)
                    {
                        if (!edr_->xdr_string(buff, MAX_STRLEN)) return TPR_FAILED;
                    }
                    break;
                }
                default:
                    THROW_TPR_EXCEPTION(
                        "Reading unknown block data type: this file is corrupted or from the "
                        "future");
            }
        }
    }

    return TPR_SUCCESS;
}

bool EdrReader::do_eheader(int nre_test)
{
    bool  bWrongPrec = false;
    bool  bOK        = false;
    float fdum       = -2e10;
    int   idum;

    // determine the precision
    const auto base = edr_->ftell_();
    if (file_version_ == 1)
    {
        edr_->fseek_(base + 12, SEEK_SET);
        int nre;
        if (!edr_->do_int(&nre)) return TPR_FAILED;
        precision_ = (nre == nre_ ? 8 : 4);
    }
    else
    {
        edr_->fseek_(base + 4, SEEK_SET);
        int magic;
        if (!edr_->do_int(&magic)) return TPR_FAILED;
        precision_ = (magic != -7777777 ? 8 : 4);
    }
    edr_->fseek_(base, SEEK_SET);
    msg("base= %lld\n", base);
    const auto dreal = (precision_ == 8 ? DataType::xdr_double : DataType::xdr_float);

    if (!edr_->do_real(&fdum, precision_)) return TPR_FAILED;
    msg("float dum= %g\n", fdum);
    if (fdum > -1e10)
    {
        // Assume we are reading an old format
        file_version_ = 1;
        // time ps
        fr_.t = fdum;
        // step
        if (!edr_->do_int(&idum)) return TPR_FAILED;
        fr_.step = idum;
    }
    else
    {
        // new edr format
        if (!edr_->do_int(&idum)) return TPR_FAILED;
        if (idum != -7777777)
        {
            THROW_TPR_EXCEPTION(
                "Energy header magic number mismatch, this is not a GROMACS edr file");
        }
        // file version
        if (!edr_->do_int(&file_version_)) return TPR_FAILED;
        if (file_version_ > edr_version)
        {
            THROW_TPR_EXCEPTION("Reading tpx file version " + std::to_string(file_version_)
                                + " with version " + std::to_string(edr_version) + " program");
        }
        // read time as double
        if (!edr_->do_double(&fr_.t)) return TPR_FAILED;
        // read step as int64
        if (!edr_->do_int64(&fr_.step)) return TPR_FAILED;
        // msg("time= %g, step= %lld\n", fr_.t, fr_.step);
        //   nsum
        if (!edr_->do_int(&fr_.nsum)) return TPR_FAILED;
        if (file_version_ >= 3)
        {
            if (!edr_->do_int64(&fr_.nsteps)) return TPR_FAILED;
        }
        else { fr_.nsteps = std::max(1, fr_.nsum); }

        if (file_version_ >= 5)
        {
            // dt as double
            if (!edr_->do_double(&fr_.dt)) return TPR_FAILED;
        }
        else { fr_.dt = 0; }
    }
    msg("fr.t= %g, fr.step= %lld\n", fr_.t, fr_.step);

    // fr.nre
    if (!edr_->do_int(&fr_.nre)) return TPR_FAILED;
    msg("fr.nre= %d\n", fr_.nre);
    // ndisre
    int ndisre = 0;
    if (file_version_ < 4)
    {
        if (!edr_->do_int(&ndisre)) return TPR_FAILED;
    }
    else
    {
        // now reserved for possible future use
        if (!edr_->do_int(&idum)) return TPR_FAILED;
    }
    msg("ndisre= %d\n", ndisre);
    // nblock
    if (!edr_->do_int(&fr_.nblock)) return TPR_FAILED;
    if (fr_.nblock < 0) { THROW_TPR_EXCEPTION("Negative nblock in edr file"); }

    if (ndisre != 0)
    {
        if (file_version_ >= 4)
        {
            THROW_TPR_EXCEPTION("Distance restraint blocks in old style in new style file");
        }
        fr_.nblock += 1;
    }

    if (nre_test >= 0
        && ((fr_.nre > 0 && fr_.nre != nre_test) || fr_.nre < 0 || ndisre < 0 || fr_.nblock < 0))
    {
        bWrongPrec = true;
        THROW_TPR_EXCEPTION("bWrongPrec");
    }

    /* we now know what these should be, or we've already bailed out because
     of wrong precision */
    if (file_version_ == 1 && (fr_.t < 0 || fr_.t > 1e20 || fr_.step < 0))
    {
        THROW_TPR_EXCEPTION(
            "edr file with negative step number or unreasonable time (and without version "
            "number).");
    }

    // msg("fr.nblock= %d\n", fr_.nblock);
    fr_.add_blocks(fr_.nblock);

    int startb = 0;
    if (ndisre > 0)
    {
        /* sub[0] is the instantaneous data, sub[1] is time averaged */
        fr_.block[0].add_subblocks(2);
        fr_.block[0].id          = enxDISRE;
        fr_.block[0].sub[0].nr   = ndisre;
        fr_.block[0].sub[1].nr   = ndisre;
        fr_.block[0].sub[0].type = dreal;
        fr_.block[0].sub[1].type = dreal;
        startb++;
    }

    for (int b = startb; b < fr_.nblock; b++)
    {
        // blocks in old version files always have 1 subblock that consists of reals.
        if (file_version_ < 4)
        {
            fr_.block[b].add_subblocks(1);
            int nrint;
            if (!edr_->do_int(&nrint)) return TPR_FAILED;
            fr_.block[b].id          = static_cast<EnumEnx>(b - startb);
            fr_.block[b].sub[0].nr   = nrint;
            fr_.block[b].sub[0].type = dreal;
        }
        else
        {
            if (!edr_->do_int(reinterpret_cast<int*>(&fr_.block[b].id))) return TPR_FAILED;
            int nsub;
            if (!edr_->do_int(&nsub)) return TPR_FAILED;
            fr_.block[b].nsub = nsub;
            fr_.block[b].add_subblocks(nsub);
            for (int i = 0; i < nsub; i++)
            {
                auto sub = &(fr_.block[b].sub[i]);
                if (!edr_->do_int(reinterpret_cast<int*>(&sub->type))) return TPR_FAILED;
                if (!edr_->do_int(&sub->nr)) return TPR_FAILED;
            }
        }
    }

    if (!edr_->do_int(&fr_.e_size)) return TPR_FAILED;

    // now reserved for possible future use
    if (!edr_->do_int(&idum)) return TPR_FAILED;
    // Do a dummy int to keep the format compatible with the old code
    if (!edr_->do_int(&idum)) return TPR_FAILED;

    if (file_version_ == 1 && nre_test < 0)
    {
        // do somethings for old edr
    }

    return TPR_SUCCESS;
}

void EdrReader::write_data(const std::string& fout) const
{
    // check data
    if (!data_.empty() && times_.size() != data_[0].value.size())
    {
        msg("%zu %zu\n", times_.size(), data_[0].value.size());
        THROW_TPR_EXCEPTION("Wrong size of times_");
    }

    FILE* fp = fopen(fout.c_str(), "w");
    if (!fp) { THROW_TPR_EXCEPTION("Can not open file to write: " + fout); }

    // write header
    fprintf(fp, "#Time/(ps),");
    for (const auto& item : data_)
    {
        fprintf(fp, "%s", item.name.c_str());
        if (!item.unit.empty()) { fprintf(fp, "/(%s)", item.unit.c_str()); }
        fprintf(fp, ",");
    }
    fprintf(fp, "\n");

    if (!data_.empty())
    {
        size_t n = data_[0].value.size();
        for (size_t i = 0; i < n; ++i)
        {
            // first write time columns
            fprintf(fp, "%g,", times_[i]);

            for (size_t j = 0; j < data_.size(); ++j)
            {
                fprintf(fp, "%g,", data_[j].value[i]);
            }
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
    fp = nullptr;
}


std::map<std::string, std::vector<double>> EdrReader::get_ene() const
{
    // check data
    if (!data_.empty() && times_.size() != data_[0].value.size())
    {
        msg("%zu %zu\n", times_.size(), data_[0].value.size());
        THROW_TPR_EXCEPTION("Wrong size of times_");
    }

    // get data
    std::map<std::string, std::vector<double>> ret;
    ret.insert(std::make_pair("Time", times_));
    for (const auto& item : data_)
    {
        ret.insert(std::make_pair(item.name, item.value));
    }

    return ret;
}
