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
      { self, ... }: {
        systems = [
          "x86_64-linux"
          "aarch64-darwin"
          "x86_64-darwin"
        ];

        perSystem =
          { pkgs, system, ... }:
          {
            packages.scripts = pkgs.stdenvNoCC.mkDerivation {
              pname = "scripts";
              version = "0.0.1";
              src = ./scripts;
              dontBuild = true;

              installPhase = ''
                mkdir -p "$out/bin"

                while IFS= read -r -d "" file; do
                  install --mode=0555 "$file" "$out/bin/"
                done < <(find . -type f -name '*.sh' -print0)
              '';
            };

            devShells.default = pkgs.mkShell {
              NIX_PATH = "nixpkgs=${pkgs.path}";

              nativeBuildInputs =
                (with pkgs; [
                  curl
                  jre
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
                      pyopenms
                    ]
                  ))
                ])
                ++ [
                  self.packages.${system}.scripts
                ];
            };
          };
      }
    );
}
