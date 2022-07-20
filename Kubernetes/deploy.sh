sudo kubectl delete job.batch/molimic-data-generation-toil --force --grace-period=0

hsdel /home/ed4bu/cath-paper-test-0-jobstore.h5
mc rb --force bournelab/home-ed4bu-cath-paper-test-0--files

if [[ $1 == "--update-docker" ]]; then
   rm -f Prop3D-0.0.1.tar
   sudo -E make -C ../Docker
   sudo docker save --output Prop3D-0.0.1.tar edraizen/Prop3D:0.0.1
   sudo k3s ctr images rm docker.io/edraizen/Prop3D:0.0.1
   sudo k3s ctr images import Prop3D-0.0.1.tar
fi

sudo kubectl apply -f deploy.yaml
sudo kubectl get all -l job-name==molimic-data-generation-toil

echo "sudo kubectl logs -f pod/$(sudo kubectl get pods -l job-name==molimic-data-generation-toil | tail -n 1 | awk '{print $1}')"
sudo kubectl logs -f pods/`sudo kubectl get pods -l job-name==molimic-data-generation-toil | tail -n 1 | awk '{print $1}'`
