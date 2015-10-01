clear all
cfg=[];
cfg.xSize=382;
cfg.ySize=390;
cfg.dataset1=1:15;
cfg.dataset2=16:30;
% cfg.scaledownup=[-10 10];
cfg.twosampletest=1; % twosampletest: 1=paired or 2=independant
cfg.backgroundfile='facebackground.tif';
cfg.rawmaps=2;
cfg.logicalmask=2 ;
% cfg.mindatapoints=15;
% cfg.nboot=1000;
imap3(cfg)


