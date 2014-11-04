import CMGTools.RootTools.fwlite.Config as cfg


tH_Yt1 = cfg.MCComponent(
    name = 'tH_Yt1',
    files = [],
    xSection = 0.01796, #PG from twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2012#MC_samples_and_cross_sections
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )


tH_YtMinus1 = cfg.MCComponent(
    name = 'tH_YtMinus1',
    files = [],
    xSection = 0.231, #PG from twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2012#MC_samples_and_cross_sections
    nGenEvents = 0,
    triggers = [],
    effCorrFactor = 1 )

tHW_Yt1 = cfg.MCComponent(
    name          = 'tHW_Yt1' ,
    files         = []        ,
    xSection      = 13*0.015*0.324*(0.215+0.0632)     , # taken from Benjamin's email
    nGenEvents    = 0         ,
    triggers      = []        ,
    effCorrFactor = 1         )

mc_tH = [
    tH_YtMinus1,
    tHW_Yt1
#    tH_Yt1,
    ]

mc_tHW = [
    tHW
    ]
