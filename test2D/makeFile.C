void makeFile()
{
    TFile * f = new TFile("/exp/icarus/app/users/nsommagg/test.root","NEW");
    TDirectory * dmu = (TDirectory*)f->mkdir("muon");
    TDirectory * dpro = (TDirectory*)f->mkdir("proton");
    TDirectory * dpion = (TDirectory*)f->mkdir("pion");
    std::array<TDirectory*,3> dirs = {dmu, dpro, dpion};
    for (auto &d : dirs)
    {
        d->cd();
        TDirectory * d1d = (TDirectory*)d->mkdir("1d");
        TDirectory * d2d = (TDirectory*)d->mkdir("2d");
        TDirectory * d2d_DNN_NOPulseTrains = (TDirectory*)d->mkdir("2d_DNN_NOPulseTrains");
        TDirectory * d2d_DNN_PulseTrains = (TDirectory*)d->mkdir("2d_DNN_PulseTrains");
        TDirectory * d1d_YZ = (TDirectory*)d->mkdir("1d_YZ");
        TDirectory * d2d_YZ = (TDirectory*)d->mkdir("2d_YZ");
        TDirectory * d2d_DNN_NOPulseTrains_YZ = (TDirectory*)d->mkdir("2d_DNN_NOPulseTrains_YZ");
        TDirectory * d2d_DNN_PulseTrains_YZ = (TDirectory*)d->mkdir("2d_DNN_PulseTrains_YZ");
    }
}