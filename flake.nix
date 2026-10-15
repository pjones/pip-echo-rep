{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    proteomics.url = "github:pjones/proteomics.nix/openms-next";
    nextflow.url = "github:pjones/nextflow.nix";
    flake-parts.url = "github:hercules-ci/flake-parts";
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } (
      { ... }: {
        systems = [
          "x86_64-linux"
          "aarch64-darwin"
          "x86_64-darwin"
        ];

        perSystem =
          { pkgs, system, ... }:
          {
            devShells.default = pkgs.mkShell {
              NIX_PATH = "nixpkgs=${pkgs.path}";

              nativeBuildInputs =
                (with pkgs; [
                  curl
                  jre
                  ruff
                  jq
                ])
                ++ (with inputs.nextflow.packages.${system}; [
                  nextflow
                  nextflow-language-server
                ])
                ++ (with inputs.proteomics.packages.${system}; [
                  comet
                  flashlfq
                  metamorpheus
                  msgfplus
                  openms
                  percolator
                  thermorawfp

                  (python3.withPackages (
                    p: with p; [
                      numpy
                      pandas
                      pipe
                      pyopenms
                      yq # For tomlq (to edit TOML files)
                    ]
                  ))
                ]);
            };
          };
      }
    );
}
