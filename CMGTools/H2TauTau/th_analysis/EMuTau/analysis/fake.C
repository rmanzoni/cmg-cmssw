void LegendSettings(TLegend *leg){
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
}

void fake(){


    gStyle->SetCanvasColor     (0);
    gStyle->SetCanvasBorderSize(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasDefH      (700);
    gStyle->SetCanvasDefW      (700);
    gStyle->SetCanvasDefX      (100);
    gStyle->SetCanvasDefY      (100);

    gStyle->SetPadColor       (0);
    gStyle->SetPadBorderSize  (10);
    gStyle->SetPadBorderMode  (0);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin   (0.08);
    gStyle->SetPadLeftMargin  (0.18);
    gStyle->SetPadRightMargin (0.1);
    gStyle->SetPadGridX       (0);
    gStyle->SetPadGridY       (0);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

    gStyle->SetLineWidth(3);
    gStyle->SetFrameFillStyle ( 0);
    gStyle->SetFrameFillColor ( 0);
    gStyle->SetFrameLineColor ( 1);
    gStyle->SetFrameLineStyle ( 0);
    gStyle->SetFrameLineWidth ( 2);
    gStyle->SetFrameBorderSize(10);
    gStyle->SetFrameBorderMode( 0);

    gStyle->SetHistFillColor(2);
    gStyle->SetHistFillStyle(0);
    gStyle->SetHistLineColor(1);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(3);

    gStyle->SetFuncColor(1);
    gStyle->SetFuncStyle(0);
    gStyle->SetFuncWidth(2);

    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(kBlack);
    gStyle->SetMarkerSize (1.4);

    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor (0);
    gStyle->SetTitleX         (0.3);

    gStyle->SetTitleSize  (0.055,"X");
    gStyle->SetTitleOffset(1.00,"X");
    gStyle->SetLabelOffset(0.005,"X");
    gStyle->SetLabelSize  (0.050,"X");
    gStyle->SetLabelFont  (42   ,"X");

    gStyle->SetLineStyleString(11,"20 10");

    gStyle->SetTitleSize  (0.055,"Y");
    gStyle->SetTitleOffset(1.400,"Y");
    gStyle->SetLabelOffset(0.010,"Y");
    gStyle->SetLabelSize  (0.050,"Y");
    gStyle->SetLabelFont  (42   ,"Y");

    gStyle->SetTextSize   (0.055);
    gStyle->SetTextFont   (42);

    gStyle->SetStatFont   (42);
    gStyle->SetTitleFont  (42);
    gStyle->SetTitleFont  (42,"X");
    gStyle->SetTitleFont  (42,"Y");

    gStyle->SetOptStat    (0);


  Double_t binning[8] = {0,10,20,30,50,70,100,200};

  Bool_t isMuon = false;
  //Bool_t isMuon = true;
  TString fname = (isMuon) ? "root_aux/Wjet_muon_training.root" : "root_aux/Wjet_electron_training.root";

  TFile *_file0 = TFile::Open(fname);

  TH1F *h_before[3];
  TH1F *h_after[3];
  TGraphAsymmErrors *gr_data[3];

  TString vname[3] = {"lepton_pt", "lepton_jetpt", "evt_njet"};
  TString vtitle[3] = {"lepton p_{T} (GeV)", "lepton jet p_{T} (GeV)", "Number of jets"};
  TCanvas *c[3];
  
  for(int ivar=0; ivar<3; ivar++){

    if(ivar==1) continue;

    TString cname = "canvas_";
    cname += vname[ivar];
    c[ivar] = new TCanvas(cname);
    
    TString hname_before = "h_bfore_";
    hname_before += vname[ivar];

    TString hname_after = "h_after_";
    hname_after += vname[ivar];

    if(ivar==0 || ivar==1){
      h_before[ivar] = new TH1F(hname_before, hname_before, 7, binning);
      h_after[ivar] = new TH1F(hname_after, hname_after, 7, binning);
    }else{
      h_before[ivar] = new TH1F(hname_before, hname_before, 10, 0, 10);
      h_after[ivar] = new TH1F(hname_after, hname_after, 10, 0, 10);
      //      h_before[ivar] = new TH1F(hname_before, hname_before, 10, -TMath::Pi(), TMath::Pi());
      //      h_after[ivar] = new TH1F(hname_after, hname_after, 10, -TMath::Pi(), TMath::Pi());
    }





    h_before[ivar]->Sumw2();
    h_after[ivar]->Sumw2();


    TString bname = vname[ivar];
    bname += " >> ";
    bname += hname_before;

    TString aname = vname[ivar];
    aname += " >> ";
    aname += hname_after;

    if(isMuon){
      kNNTrainingTree->Draw(bname,"(evt_nbjet>=1 && (!evt_isMC || evt_id==0 || evt_id==1 || evt_id==24 || evt_id==25))*evt_weight*evt_isMCw");
      kNNTrainingTree->Draw(aname, "(evt_nbjet>=1 && (!evt_isMC || evt_id==0 || evt_id==1 || evt_id==24 || evt_id==25) && lepton_id==1 && lepton_mva > lepton_mva_threshold)*evt_weight*evt_isMCw");
    }else{
      kNNTrainingTree->Draw(bname,"(evt_nbjet>=0 && (!evt_isMC || evt_id==0 || evt_id==1 || evt_id==24 || evt_id==25))*evt_weight*evt_isMCw");
      kNNTrainingTree->Draw(aname, "(evt_nbjet>=0 && (!evt_isMC || evt_id==0 || evt_id==1 || evt_id==24 || evt_id==25) && lepton_id==1 && lepton_mva > lepton_mva_threshold)*evt_weight*evt_isMCw");
    }
    
    gr_data[ivar] = new TGraphAsymmErrors();

    gr_data[ivar]->BayesDivide(h_after[ivar],h_before[ivar]);   
    gr_data[ivar]->SetMarkerStyle(20);
    gr_data[ivar]->SetMarkerSize(1);
    //    gr_data[ivar]->GetXaxis()->SetTitle("lepton p_{T} (GeV)");
    gr_data[ivar]->GetXaxis()->SetTitle(vtitle[ivar]);
    gr_data[ivar]->GetYaxis()->SetTitle("Fake Rate");
    gr_data[ivar]->GetXaxis()->SetNdivisions(505);
    gr_data[ivar]->SetMarkerStyle(20);
    gr_data[ivar]->SetLineWidth(2);
    gr_data[ivar]->SetMaximum(0.3);
    gr_data[ivar]->SetMinimum(0.);
    gr_data[ivar]->Draw("aep");

    TString lname = (isMuon) ? "jet #rightarrow #mu" : "jet #rightarrow e";

    Float_t xval = (ivar==0) ? 152 : 7.5;
    TLatex * tex = new TLatex(xval, 0.2644593,lname);
    tex->SetTextFont(42);
    tex->SetTextSize(0.05637982);
    tex->SetLineWidth(2);
    tex->Draw();

    TString savename = vname[ivar];
    savename += (isMuon) ? "_muon" : "_electron";
    savename += ".pdf";
    c[ivar]->SaveAs(savename);

    TString savename2 = vname[ivar];
    savename2 += (isMuon) ? "_muon" : "_electron";
    savename2 += ".gif";
    c[ivar]->SaveAs(savename2);
  }

  

    



}
