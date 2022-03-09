//
// Fast Cross-Correlation algorithm
//
// Author: Francois Grondin
// Email: francois.grondin2@usherbrooke.ca
//
// Refer to the following paper for details:
//
// Grondin, F., Maheux, M.-A., Lauzon, J.-S. and Michaud, F., "Fast Cross-Correlation
// for TDoA Estimation on Small Aperture Microphone Arrays", arXiV
//

#include <fccphat/system.h>

#include <gtest/gtest.h>

#define EXPECT_EQ_ARRAY(val1, val2, size) do {         \
        for (size_t i = 0; i < size; i++) {         \
            EXPECT_EQ(val1[i], val2[i]) << "(index=" << i << ")"; \
        }                                           \
    } while(false)

#define EXPECT_NEAR_ARRAY(val1, val2, size, abs_error) do { \
        for (size_t i = 0; i < size; i++) {              \
            EXPECT_NEAR(val1[i], val2[i], abs_error) << "(index=" << i << ")";    \
        }                                                \
    } while(false)

TEST(system_tests, wav_obj_should_be_able_to_read_a_wav_file) {
    wav_obj * wav = nullptr;
    hops_obj * hops = nullptr;

    const char * wav_filename = "resources/test_wav_obj.wav";
    const unsigned int hop_size = 3;

    const int16_t expected_buffer_1[] = {327, 2293, 655, 2621, 983, 2949};
    const int16_t expected_buffer_2[] = {1310, 3276, 1638, 3604, 1966, 3932};
    const float expected_samples_1[] = {0.01f, 0.02f, 0.03f};
    const float expected_samples_2[] = {0.07f, 0.08f, 0.09f};
    const float expected_samples_3[] = {0.04, 0.05f, 0.06f};
    const float expected_samples_4[] = {0.1f, 0.11f, 0.12f};
    const float samples_abs_error = 0.001;


    wav = wav_construct(wav_filename, hop_size);

    EXPECT_NE(wav, nullptr);
    EXPECT_EQ(wav->sample_rate, 44100);
    EXPECT_EQ(wav->num_channels, 2);
    EXPECT_EQ(wav->bits_per_sample, 16);
    EXPECT_EQ(wav->hop_size, 3);

    hops = hops_construct(wav->num_channels, hop_size);
    EXPECT_EQ(hops->channels_count, 2);
    EXPECT_EQ(hops->hop_size, 3);
    EXPECT_NE(hops, nullptr);

    wav_read(wav, hops);
    EXPECT_EQ_ARRAY(wav->buffer, expected_buffer_1, wav->num_channels * wav->hop_size);
    EXPECT_NEAR_ARRAY(hops->samples[0], expected_samples_1, wav->hop_size, samples_abs_error);
    EXPECT_NEAR_ARRAY(hops->samples[1], expected_samples_2, wav->hop_size, samples_abs_error);

    wav_read(wav, hops);
    EXPECT_EQ_ARRAY(wav->buffer, expected_buffer_2, wav->num_channels * wav->hop_size);
    EXPECT_NEAR_ARRAY(hops->samples[0], expected_samples_3, wav->hop_size, samples_abs_error);
    EXPECT_NEAR_ARRAY(hops->samples[1], expected_samples_4, wav->hop_size, samples_abs_error);

    hops_destroy(hops);
    wav_destroy(wav);
}
