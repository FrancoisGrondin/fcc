#include <fccphat/signal.h>
#include <fccphat/system.h>

#include <unistd.h>
#include <string.h>

int main(int argc, char * argv[]) {

    int opt;
    char method;
    char filename[1024];

    wav_obj * wav;
    stft_obj * stft;
    scmphat_obj * scmphat;
    gcc_obj * gcc;
    fcc_obj * fcc;
    quadinterp_obj * quadinterp;

    hops_obj * hops;
    freqs_obj * freqs;
    covs_obj * covs;
    corrs_obj * corrs;
    taus_obj * taus;

    const unsigned int hop_size = 256;
    const unsigned int frame_size = 512;
    const unsigned int sample_rate = 16000;
    const unsigned int tau_max = 8;
    const unsigned int interpolation_rate = 2;
    const float alpha = 1.0;

    method = 0x00;

    while((opt = getopt(argc, argv, "i:m:")) != -1) 
    { 
        switch(opt) 
        { 
            
            case 'i':

                strcpy(filename, optarg);
                break;

            case 'm':

                if (strcmp(optarg, "gcc") == 0) { method = 'g'; }
                if (strcmp(optarg, "fcc") == 0) { method = 'f'; }
                break;

        } 
    } 

    if (method == 0x00) {
        printf("Invalid method. Choose either \"fcc\" or \"gcc\"\n");
        exit(-1);        
    }

    if( access( filename, F_OK ) != 0 ) {
        printf("File %s does not exists\n", filename);
        exit(-1);
    }

    wav = wav_construct(filename, hop_size);

    if (wav->sample_rate != sample_rate) {
        printf("Invalid sample rate (was expecting %u, got %u).\n", sample_rate, wav->sample_rate);
        exit(-1);
    }


    stft = stft_construct(wav->num_channels, frame_size, hop_size);
    scmphat = scmphat_construct(wav->num_channels, frame_size, alpha);
    gcc = gcc_construct(wav->num_channels, frame_size, tau_max, interpolation_rate);
    fcc = fcc_construct(wav->num_channels, frame_size, tau_max);
    quadinterp = quadinterp_construct(wav->num_channels);

    hops = hops_construct(wav->num_channels, hop_size);
    freqs = freqs_construct(wav->num_channels, frame_size);
    covs = covs_construct(wav->num_channels, frame_size);
    corrs = corrs_construct(wav->num_channels);
    taus = taus_construct(wav->num_channels);

    while (wav_read(wav, hops) == 0) {

        stft_call(stft, hops, freqs);
        scmphat_call(scmphat, freqs, covs);

        if (method == 'g') {
            gcc_call(gcc, covs, corrs);
        }
        if (method == 'f') {
            fcc_call(fcc, covs, corrs);
        }
        
        quadinterp_call(quadinterp, corrs, taus);

    }

    stft_destroy(stft);
    scmphat_destroy(scmphat);
    gcc_destroy(gcc);
    fcc_destroy(fcc);
    quadinterp_destroy(quadinterp);

    hops_destroy(hops);
    freqs_destroy(freqs);
    covs_destroy(covs);
    corrs_destroy(corrs);
    taus_destroy(taus);

    return 0;

}