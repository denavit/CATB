function Liu_et_al_2013_examples()
% Liu et al. 2013 - Examples

fid = 1;
% fid = fopen('Liu_et_al_2013_examples_results.txt','w');

%% Example 1
wf = wf_caftb('W14x48',50,'AISC2010');
fprintf(fid,'\nExample 1\n');
print_output(fid,368,wf.phi_c*wf.Pnz(32*12,1));

%% Example 2
wf = wf_caftb('W16x26',50,'AISC2010');
fprintf(fid,'\nExample 2\n');
print_output(fid,235,wf.phi_c*wf.Pnz(8*12,1));

%% Example 3
wf = wf_caftb('W14x90',50,'AISC2010');
fprintf(fid,'\nExample 3\n');
print_output(fid,1000,wf.phi_c*wf.Pny(15*12,1));
print_output(fid, 928,wf.phi_c*wf.Pnx(30*12,1));
print_output(fid, 838,wf.phi_c*wf.Pnz(30*12,1));

%% Example 4
wf = wf_caftb('W18x35',50,'AISC2010');
fprintf(fid,'\nExample 4\n');
print_output(fid,299,wf.phi_c*wf.Pnca(8*12,1));

%% Example 5
wf = wf_caftb('W14x132',50,'AISC2010');
fprintf(fid,'\nExample 5\n');
print_output(fid,1090,wf.phi_c*wf.Pnca(40*12,1));

%% 
if fid ~= 1
    fclose(fid)
end

end

function print_output(fid,paper_answer,calculated_answer)
    fprintf(fid,'%.0f (Paper) vs. %.1f (Calculated)\n',paper_answer,calculated_answer);
    percent_diff = 100*(paper_answer-calculated_answer)/paper_answer;
    if abs(percent_diff) < 1
        fprintf(fid,'   OKAY (Percent Difference = %.1f%%)\n',percent_diff);
    else
        fprintf(fid,'   CHECK RESULT (Percent Difference = %.1f%%)\n',percent_diff);
    end
end