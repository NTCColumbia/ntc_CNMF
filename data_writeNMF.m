function [ output_args ] = data_writeNMF(Datawrite,weigh_image,A, Ain, C, Cin, b, f, signal_raw, signal_filtered, signal_inferred,Y_fres,use_merged  )
%data_write_NMF(Datawrite,A, Ain, C, Cin, b, f, signal_raw, signal_filtered, signal_inferred)
%Horribly kludged quick function to allow saving of cNMF output parameters
%and traces into a unique file name and structure to allow for "easy"
%reloading into matlab and file comparisons.  Many better ways to deal with
%this, but this is what I have today...6/10/15 Darcy S. Peterka

appe=datetime;
appe=strrep(char(appe), ' ', '_');
appe=strrep(char(appe), ':', '_');
appe=strrep(char(appe), '-', '_');
tmpstruct.datawrite=Datawrite;
tmpstruct.A=A;
tmpstruct.Ain=Ain;
tmpstruct.C=C;
tmpstruct.Cin=Cin;
tmpstruct.b=b;
tmpstruct.f=f;
tmpstruct.signal_raw=signal_raw;
tmpstruct.signal_filtered=signal_filtered;
tmpstruct.signal_inferred=signal_inferred;
tmpstruct.weigh_image=weigh_image;
tmpstruct.Y_fres=Y_fres;

[~,tmpjunkn, ~]=fileparts(tmpstruct.datawrite.movieFileName);
if use_merged==0;
tmpjnk=strcat(tmpjunkn,'_merged_')
savedstruct=strcat(tmpjnk,appe);
savedfile=strcat(tmpjnk,'_merged_',appe,'_cNMF.mat');
else
savedstruct=strcat(tmpjunkn,appe);
savedfile=strcat(tmpjunkn,appe,'_cNMF.mat');
end
vtmp=genvarname(savedstruct);
eval([vtmp '= tmpstruct;']);
save(savedfile,savedstruct);
datawritestr=['Wrote out file ', [savedfile], ' to disk'];
    display(datawritestr);
clear tmpstruct appe savedstruct savedfile vtmp
end

