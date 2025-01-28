function [Y, H] = get_data_for_detector(m, k, beam_amplitudes, h_data_recovered_f, h_data_noisy, transform_matrix, scenario)

if ~scenario.half_precision
    switch scenario.beam_transform
        case 1
            H=squeeze(beam_amplitudes(:,m,:,k));
            if scenario.N_target_UE==1
                H=H.';
            end
            Y=(squeeze(h_data_noisy(:,m,k))).'*conj(transform_matrix);
        case 0
            H=squeeze(h_data_recovered_f(:,:,m,k));
            Y=squeeze(h_data_noisy(:,m,k)).';
    end
    Y = single(Y.');
    H = single(H.');
elseif scenario.half_precision

    sc_range = m;
    N_used = length(sc_range);

    switch scenario.beam_transform
        case 1
            H=permute(beam_amplitudes(:,sc_range,:,:), [3,1,4,2]);
%             if scenario.N_target_UE==1
%                 H=reshape(H, [1,scenario.N_ports,scenario.N_data_sym]);
%                 %H=permute(H, [3,1,2]);
%             end
            Y = zeros(scenario.N_data_sym, N_used, scenario.N_ports);
            for sc = 1:length(sc_range)
                Y(:,sc,:)=(squeeze(h_data_noisy(:,sc_range(sc),:))).'*conj(transform_matrix);
            end
            Y = permute(Y, [3,1,2]);
        case 0
            H=squeeze(h_data_recovered_f(:,:,m,:));
            if scenario.N_target_UE==1
                H=reshape(H, [1,scenario.Nrx,scenario.N_data_sym, ]);
            end
            Y= permute(h_data_noisy, [1,3,2]);
    end

    Y = single(Y);
    H = single(H);
end