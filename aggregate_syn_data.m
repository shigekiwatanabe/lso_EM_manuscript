%ALL VARIABLES
    %syn_vesicles - structure of vesicle measurements from each bouton (each row is one synapse)
        %count - number of vesicles per bouton
        %major - array containing major axes of vesicles 
        %minor - array containing minor axes of vesicles 
        %area - array containing areas of all vesicles 
        %ar - array containing aspect ratios of vesicles
        %stde - standard deviation of vesicle aspect ratio
        %cv - coefficient of variation of vesicle aspect ratio
        %ar_avg - average aspect ratio
        %exc_count - number of vesicles (excitatory only)
        %exc_area - array containing areas of vesicles (excitatory only)
        %inh_count - number of vesicles (inhibitory only)
        %inh_area - array containing areas of vesicles (inhibitory only)
    %mito - structure with mitochondrion measurements from each bouton (each row is one bouton)
        %count - number of mitochondrion per bouton
        %area - array containing areas of all mitochondria
    %syn - structure of bouton measurements (each row is one bouton)
        %area - area of bouton
        %density - vesicle density
        %contact - length of bouton membrane contact with postsynaptic cell
        %az - array containing length of active zones
        %total_az - summed length of all active zones 
    %exc_syn - structure of excitatory bouton measurements
        %area, density, az
    %inh_syn - structure of inhibitory bouton measurements
        %area, density, az
    %nn - array containing distance of 3 nearest neighboring vesicles from each vesicle   
    %inh_nn - nearest neighbors for inhibitory boutons
    %exc_nn - nearest neighbors for excitatory boutons
    %file - bouton file name

function y = aggregate_syn_data (x)
    inh_count = 0;
    exc_count = 0;

    for i = 1:length(x.raw_data)
        %vesicles
        if isfield(x.raw_data(i).analysis_data, "vesicle") == 1
            num_vesicles = length(x.raw_data(i).analysis_data.vesicle);
            y.syn_vesicles(i).count = num_vesicles;

            for k = 1:num_vesicles
                y.syn_vesicles(i).major(k) = x.raw_data(i).analysis_data.vesicle(k).major;
                y.syn_vesicles(i).minor(k) = x.raw_data(i).analysis_data.vesicle(k).minor;
                y.syn_vesicles(i).area(k) = x.raw_data(i).analysis_data.vesicle(k).areapx;
                y.syn_vesicles(i).ar(k) = x.raw_data(i).analysis_data.vesicle(k).ar;
            end

            %analyze average AR for synapses with more than 5 vesicles
            if num_vesicles >=5
                y.syn_vesicles(i).stde = std(y.syn_vesicles(i).ar,0,'all');
                y.syn_vesicles(i).cv = std(y.syn_vesicles(i).ar,0,'all')/mean(y.syn_vesicles(i).ar,'all');
                y.syn_vesicles(i).ar_avg = mean(y.syn_vesicles(i).ar);
            else
                y.syn_vesicles(i).exc_count = [];
                y.syn_vesicles(i).inh_count = [];
                y.syn_vesicles(i).ar_avg = [];
            end

            %vesicle parameters for inhibitory and excitatory synapses
            %(identified by threshold of 1.33 AR)
            if (y.syn_vesicles(i).ar_avg >= 1.33)
                inh_count = inh_count + 1;
                y.syn_vesicles(i).inh_count = num_vesicles;
                y.syn_vesicles(i).inh_area = y.syn_vesicles(i).area;
    
                y.syn_vesicles(i).exc_count = [];
                y.syn_vesicles(i).exc_area = [];
            end
    
            if (y.syn_vesicles(i).ar_avg < 1.33)
                exc_count = exc_count + 1;
                y.syn_vesicles(i).exc_count = num_vesicles;
                y.syn_vesicles(i).exc_area = y.syn_vesicles(i).area;
    
                y.syn_vesicles(i).inh_count = [];
                y.syn_vesicles(i).inh_area = [];
            end

        else
            %no vesicles
            y.syn_vesicles(i).count = 0;
            y.syn_vesicles(i).major = [];
            y.syn_vesicles(i).minor = [];
            y.syn_vesicles(i).area = [];
            y.syn_vesicles(i).cv = [];
            y.syn_vesicles(i).ar = [];
    
        end

        %mitochondria
        if isfield(x.raw_data(i).analysis_data, 'endosome') == 1
            y.mito(i).count = length(x.raw_data(i).analysis_data.endosome);
            for j = 1:y.mito(i).count
                y.mito(i).area(j) = x.raw_data(i).analysis_data.endosome(j).area;
            end   
        else
            y.mito(i).count = 0;
            y.mito(i).area = NaN;
        end
    
        %plasma membrane and density
        if isfield(x.raw_data(i).analysis_data.pm, 'area') == 1
            y.syn(i).area = x.raw_data(i).analysis_data.pm.area;
            y.syn(i).density = y.syn_vesicles(i).count/y.syn(i).area;
            if y.syn_vesicles(i).count >= 3
                y.syn(i).density3 = y.syn(i).density;
            end
            if isfield(x.raw_data(i).analysis_data, 'az') == 1
                y.syn(i).contact = x.raw_data(i).analysis_data.az.length;
            else
                y.syn(i).contact = [];
            end
            if isfield(x.raw_data(i).analysis_data,'dp') == 1
                y.syn(i).az = [x.raw_data(i).analysis_data.dp.length];
                y.syn(i).total_az = sum([y.syn(i).az]);
            else
                y.syn(i).az = [];
            end
    
            if isfield(x.raw_data(i).analysis_data, "vesicle") == 1
                if (y.syn_vesicles(i).ar_avg >= 1.33)
                    y.exc_syn(i).area = [];
                    y.inh_syn(i).area = y.syn(i).area;
                    y.exc_syn(i).density = [];
                    y.inh_syn(i).density = y.syn(i).density;
                    y.exc_syn(i).az = [];
                    y.inh_syn(i).az = y.syn(i).az;
                elseif (y.syn_vesicles(i).ar_avg < 1.33)
                    y.exc_syn(i).area = y.syn(i).area;
                    y.inh_syn(i).area = [];
                    y.exc_syn(i).density = y.syn(i).density;
                    y.inh_syn(i).density = [];
                    y.exc_syn(i).az = y.syn(i).az;
                    y.inh_syn(i).az = [];
                else
                    y.inh_syn(i).area = [];
                    y.inh_syn(i).density = [];
                    y.inh_syn(i).az = [];
                    y.exc_syn(i).area = [];
                    y.exc_syn(i).density = [];
                    y.exc_syn(i).az = [];
                end
            end
        end
    %save file name
    y.file(i).name = x.raw_data(i).file_name;
    end

    %aggregate vesicle parameters
    c = 1;
    inh_v_count = 1;
    exc_v_count = 1;
    for i = 1: length(y.syn_vesicles)
        for j = 1: length(y.syn_vesicles(i).area)
            y.all_sv(c).vesicles_area = y.syn_vesicles(i).area(j);
            y.all_sv(c).vesicles_major = y.syn_vesicles(i).major(j);
            y.all_sv(c).vesicles_minor = y.syn_vesicles(i).minor(j);
            y.all_sv(c).vesicles_ar = y.syn_vesicles(i).ar(j);
            c = c + 1;
            if (y.syn_vesicles(i).ar_avg >= 1.33)
                y.all_sv(inh_v_count).inh_vesicles_area = y.syn_vesicles(i).area(j);
                y.all_sv(inh_v_count).inh_vesicles_major = y.syn_vesicles(i).major(j);
                y.all_sv(inh_v_count).inh_vesicles_minor = y.syn_vesicles(i).minor(j);
                y.all_sv(inh_v_count).inh_vesicles_ar = y.syn_vesicles(i).ar(j);
                inh_v_count = inh_v_count+ 1;
            end
            if (y.syn_vesicles(i).ar_avg < 1.33)
                y.all_sv(exc_v_count).exc_vesicles_area = y.syn_vesicles(i).area(j);
                y.all_sv(exc_v_count).exc_vesicles_major = y.syn_vesicles(i).major(j);
                y.all_sv(exc_v_count).exc_vesicles_minor = y.syn_vesicles(i).minor(j);
                y.all_sv(exc_v_count).exc_vesicles_ar = y.syn_vesicles(i).ar(j);
                exc_v_count = exc_v_count+ 1;
            end
        end
    end


    %excitatory and inhibitory mitochondria
    inh_m_count = 1;
    exc_m_count = 1;
    m = 1;
    
    for i = 1:length(y.mito)
        for j = 1:length(y.mito(i).area)
            if y.mito(i).count ~= 0 && isfield(x.raw_data(i).analysis_data, "vesicle") == 1
            y.all_mito(m).area = y.mito(i).area(j);
                if (y.syn_vesicles(i).ar_avg >= 1.33)
                    y.all_mito(inh_m_count).inh_mito_area = y.all_mito(m).area;
                    inh_m_count = inh_m_count + 1;
                elseif (y.syn_vesicles(i).ar_avg < 1.33)
                    y.all_mito(exc_m_count).exc_mito_area = y.all_mito(m).area;
                    exc_m_count = exc_m_count + 1;
                end
            m = m+1;
            end
            
        end
    end

%nearest neighbors
[y.nn, y.inh_nn, y.exc_nn] = find_nn3 (x,y);

end

