function [all_nn, all_inh_nn, all_exc_nn] = find_nn3 (x,y)
    for i = 1:length(x.raw_data)
        if isfield(x.raw_data(i).analysis_data, "vesicle") == 1 
            %vesicle x, y input
            v = [x.raw_data(i).analysis_data.vesicle.x; x.raw_data(i).analysis_data.vesicle.y];
            v_ar = [x.raw_data(i).analysis_data.vesicle.ar];
            
    
            %transpose
            vT = v.';
            vT_ar = v_ar;
            
            if height(vT) > 4
                for j = 1:length(vT)
                    %search array
                    s = vT;
                    s_ar = vT_ar;
                    s(j,:) = [];
                    s_ar(j) = [];
                    
                    %query array
                    q = vT(j,:);
                    
                    %nearest neighbor
                    [Idx, D] = knnsearch(s, q, 'Distance', 'euclidean', 'K', 3);
                    %nn(i).idx = Idx;
                    nn(i).v(j).dist = D;
                    nn(i).v(j).ar = s_ar(Idx);
                    nn(i).v(j).ar_diff = vT_ar(j) - s_ar(Idx);
                end
            end
        end
    end
    
    %all nearest neighbors
    c = 1;
    inh_c = 1;
    exc_c = 1;
    for i = 1:length(nn)
        if length(nn(i).v) ~= 0
            for j = 1:length(nn(i).v)
                if isfield(nn(i).v(j), 'dist') == 1
                    all_nn1(c) = nn(i).v(j).dist(1);
                    all_nn2(c) = nn(i).v(j).dist(2);
                    all_nn3(c) = nn(i).v(j).dist(3);
                    c = c+1;
                    if y.syn_vesicles(i).ar_avg >= 1.33
                        all_inh_nn1(inh_c) = nn(i).v(j).dist(1);
                        all_inh_nn2(inh_c) = nn(i).v(j).dist(2);
                        all_inh_nn3(inh_c) = nn(i).v(j).dist(3);
                        inh_c = inh_c + 1;
                    end
                    if y.syn_vesicles(i).ar_avg < 1.33
                        all_exc_nn1(exc_c) = nn(i).v(j).dist(1);
                        all_exc_nn2(exc_c) = nn(i).v(j).dist(2);
                        all_exc_nn3(exc_c) = nn(i).v(j).dist(3);
                        exc_c = exc_c + 1;
                    end
                end
                
            end
        end
    end
    
    all_nn = [all_nn1', all_nn2', all_nn3'];
    all_inh_nn = [all_inh_nn1', all_inh_nn2', all_inh_nn3'];
    all_exc_nn = [all_exc_nn1', all_exc_nn2', all_exc_nn3'];

end
%separate by exc inh

