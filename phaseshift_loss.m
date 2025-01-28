function loss = phaseshift_loss(exps, H, inds, uniformity)
    N_digital = size(inds, 2);
    N_rx = size(H,1);

    Fanalog = zeros(N_rx, N_digital);
    for i = 1:N_digital
        Fanalog(inds(:,i), i) = exps(inds(:,i));
    end
    Projection = Fanalog'*H;
    power_per_sc = sum(abs(Projection).^2,1);
%     power_per_sc = sum(abs(Projection),1);
%     warning("No second power of projection!")
    
    loss = - sum(power_per_sc);
    if uniformity
        term_uniformity = sum(abs(power_per_sc - mean(power_per_sc)));
        loss = loss + term_uniformity;
    end
end
    