import React from 'react';
import ReactCollapsible from 'react-collapsible';
import './Collapsible.style.scss';

const CollapsibleTrigger = (props) => {
  const { trigger } = props;
  return (
    <div className='collapsible-trigger--inner'>
      {trigger}
      <span className='material-icons trigger-icon'>keyboard_arrow_down</span>
    </div>
  );
};

const Collapsible = (props) => {
  const { trigger, children } = props;
  const reactCollapsibleRef = React.createRef<ReactCollapsible>();
  const handleCollapsibleOpen = () => {
    console.log(
      ((reactCollapsibleRef.current as any).innerRef.style.overflow = 'visible')
    );
  };
  const handleCollapsibleClosing = () => {
    (reactCollapsibleRef.current as any).innerRef.style.overflow = 'hidden';
  };
  return (
    <div className='collapsible--custom'>
      <ReactCollapsible
        trigger={<CollapsibleTrigger trigger={trigger} />}
        onOpen={handleCollapsibleOpen}
        onClosing={handleCollapsibleClosing}
        ref={reactCollapsibleRef}
      >
        {children}
      </ReactCollapsible>
    </div>
  );
};

export default Collapsible;
